"""Write the small test files in tests/fixtures and their expected slices.

Requires numpy and h5py. Expected values come from average_slab(), copied
unchanged from the original Python slice viewer, so the JavaScript port is
checked against the reference implementation.
"""
import json
from pathlib import Path

import h5py
import numpy as np

OUT = Path(__file__).resolve().parents[1] / 'tests' / 'fixtures'


def average_slab(group, fixed, indices, inversion=False):
    """Pool unique finite source voxels across the slab and inversion partner.

    Equal weight per source voxel. Errors assume independent source voxels;
    covariance introduced by prior data reduction is not available here.
    Returned arrays follow the remaining storage axes (storage order L,K,H).
    """
    shape = group["signal"].shape
    axis = 2 - fixed
    remaining = [d for d in range(3) if d != axis]
    grids = np.meshgrid(*(np.arange(shape[d]) for d in remaining), indexing="ij")
    total = np.zeros(grids[0].shape)
    variance_sum = np.zeros_like(total)
    count = np.zeros(total.shape, dtype=np.int16)
    error_complete = np.ones(total.shape, dtype=bool)
    seen = []
    for inverted in ([False, True] if inversion else [False]):
        for index in indices:
            source_index = shape[axis] - 1 - index if inverted else index
            sel = [slice(None)] * 3
            sel[axis] = int(source_index)
            values = np.asarray(group["signal"][tuple(sel)])
            mask = np.asarray(group["mask"][tuple(sel)])
            variances = np.asarray(group["errors_squared"][tuple(sel)])
            if inverted:
                values, mask, variances = (v[::-1, ::-1] for v in (values, mask, variances))
            coords = [None] * 3
            coords[axis] = np.full(total.shape, source_index, dtype=int)
            for d, grid in zip(remaining, grids):
                coords[d] = shape[d] - 1 - grid if inverted else grid
            source_ids = np.ravel_multi_index(tuple(coords), shape)
            unique = np.ones(total.shape, dtype=bool)
            for previous in seen:
                unique &= source_ids != previous
            seen.append(source_ids)
            valid = unique & np.isfinite(values) & (mask == 0)
            known_error = np.isfinite(variances) & (variances >= 0)
            total += np.where(valid, values, 0)
            variance_sum += np.where(valid & known_error, variances, 0)
            error_complete &= ~valid | known_error
            count += valid
    mean = np.divide(total, count, out=np.full(total.shape, np.nan), where=count > 0)
    variance = np.divide(variance_sum, count.astype(float)**2,
                         out=np.full(total.shape, np.nan), where=(count > 0) & error_complete)
    return mean, count, variance


def orientation_matrix(a, b, c, alpha, beta, gamma):
    """UB = U B with B^T B the reciprocal metric (no 2*pi) and a fixed rotation U."""
    ca, cb, cg = np.cos(np.radians([alpha, beta, gamma]))
    g = np.array([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]])
    B = np.linalg.cholesky(np.linalg.inv(g)).T
    t = np.radians(35)
    U = np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])
    return U @ B


def cases(signal, mask, edges, specs):
    group = dict(signal=signal, mask=mask, errors_squared=np.broadcast_to(0., signal.shape))
    centers = [(e[1:] + e[:-1]) / 2 for e in edges]
    out = []
    for fixed, center, thickness, inversion in specs:
        ids = np.flatnonzero(abs(centers[fixed] - center) < thickness / 2 + 1e-5)
        z, count, _ = average_slab(group, fixed, ids, inversion)
        out.append(dict(fixed=fixed, center=center, thickness=thickness, inversion=inversion,
                        ids=ids.tolist(), shape=list(z.shape),
                        values=[None if np.isnan(v) else float(v) for v in z.ravel()],
                        counts=count.ravel().tolist()))
    return out


def mdhisto(rng):
    """Mantid SaveMD layout: storage (L, K, H), origin-centered edges, hexagonal cell."""
    shape = (7, 9, 11)
    signal = rng.normal(10, 4, shape)
    signal[rng.random(shape) < 0.08] = np.nan
    signal[rng.random(shape) < 0.03] = np.inf
    mask = (rng.random(shape) < 0.12).astype(np.int8)
    edges = [np.linspace(-(n / 2) * w, (n / 2) * w, n + 1) for n, w in zip(shape[::-1], (0.1, 0.1, 0.25))]
    with h5py.File(OUT / 'mdhisto_small.nxs', 'w') as f:
        entry = f.create_group('MDHistoWorkspace')
        entry.attrs['NX_class'] = 'NXentry'
        data = entry.create_group('data')
        data.attrs['NX_class'] = 'NXdata'
        kw = dict(chunks=(1,) + shape[1:], compression='gzip', compression_opts=4)
        ds = data.create_dataset('signal', data=signal, **kw)
        ds.attrs['signal'] = 1
        ds.attrs['axes'] = np.bytes_('D2:D1:D0')
        data.create_dataset('mask', data=mask, **kw)
        data.create_dataset('errors_squared', data=np.abs(signal), **kw)
        for i, (name, e) in enumerate(zip(['[H,0,0]', '[0,K,0]', '[0,0,L]'], edges)):
            d = data.create_dataset(f'D{i}', data=e)
            d.attrs.update(frame=np.bytes_('HKL'), long_name=np.bytes_(name), units=np.bytes_('r.l.u.'))
        cell = entry.create_group('experiment0/sample/oriented_lattice')
        cell.attrs['NX_class'] = 'NXcrystal'
        params = dict(a=8.03, b=8.02, c=10.03, alpha=90.1, beta=90.2, gamma=119.9)
        for k, v in params.items():
            cell.create_dataset(f'unit_cell_{k}', data=[v])
        cell.create_dataset('orientation_matrix', data=orientation_matrix(**params))
    specs = [(2, 0.0, 0.3, True), (2, 0.0, 0.3, False), (1, 0.1, 0.25, True),
             (0, 0.0, 0.2, True), (0, -0.4, 0.5, False), (2, 0.25, 1.0, True), (1, 0.0, 10.0, True)]
    return cases(signal, mask, edges, specs)


def nxdata(rng):
    """Plain NXdata: 4-D with a size-1 axis, bin-center axes, float32, no mask or cell."""
    shape = (1, 6, 8, 10)
    signal = rng.gamma(2, 3, shape).astype(np.float32)
    signal[rng.random(shape) < 0.1] = np.nan
    axes = dict(t=np.array([5.0]), qz=np.linspace(0.5, 3.0, 6), qy=np.linspace(-1.4, 1.4, 8),
                qx=np.linspace(-0.9, 0.9, 10))
    with h5py.File(OUT / 'nxdata_small.h5', 'w') as f:
        entry = f.create_group('entry')
        entry.attrs['NX_class'] = 'NXentry'
        data = entry.create_group('data')
        data.attrs.update(NX_class='NXdata', signal='intensity', axes=['t', 'qz', 'qy', 'qx'])
        data.create_dataset('intensity', data=signal)
        for name, values in axes.items():
            data.create_dataset(name, data=values).attrs['units'] = 'Angstrom^-1'
    squeezed = signal[0].astype(np.float64)

    def to_edges(c):
        mid = (c[1:] + c[:-1]) / 2
        return np.concatenate([[c[0] - (c[1] - c[0]) / 2], mid, [c[-1] + (c[-1] - c[-2]) / 2]])
    edges = [to_edges(axes[k]) for k in ('qx', 'qy', 'qz')]
    specs = [(2, 1.5, 0.4, False), (1, 0.2, 0.5, False), (0, 0.0, 0.3, False)]
    return cases(squeezed, np.zeros(squeezed.shape, np.int8), edges, specs)


def centers_to_edges(c):
    mid = (c[1:] + c[:-1]) / 2
    return np.concatenate([[c[0] - (c[1] - c[0]) / 2], mid, [c[-1] + (c[-1] - c[-2]) / 2]])


def nebula(rng):
    """NEBULA3D volume: /entry, storage (H, K, L), bin centers, mask True = valid, sigma, UB with 2*pi."""
    shape = (5, 7, 9)
    signal = rng.normal(5, 2, shape)
    signal[rng.random(shape) < 0.05] = np.nan
    valid = rng.random(shape) > 0.1
    sigma = np.abs(rng.normal(0.5, 0.1, shape))
    axes = dict(h_axis=np.linspace(-0.4, 0.4, 5), k_axis=np.linspace(-0.3, 0.3, 7), l_axis=np.linspace(-1, 1, 9))
    params = dict(a=5.91, b=10.42, c=24.79, alpha=89.55, beta=90.61, gamma=90.63)
    with h5py.File(OUT / 'nebula3d_small.h5', 'w') as f:
        entry = f.create_group('entry')
        kw = dict(compression='gzip', compression_opts=1, shuffle=True)
        entry.create_dataset('data', data=signal, **kw)
        entry.create_dataset('sigma', data=sigma, **kw)
        entry.create_dataset('mask', data=valid, **kw)
        for name, values in axes.items():
            entry.create_dataset(name, data=values)
        entry.create_dataset('ub_matrix', data=2 * np.pi * orientation_matrix(**params))
        entry.attrs['instrument'] = 'synthetic'
    # Displayed fastest storage axis first: L, K, H.
    edges = [centers_to_edges(axes[k]) for k in ('l_axis', 'k_axis', 'h_axis')]
    specs = [(2, 0.0, 0.3, True), (1, 0.1, 0.25, False), (0, -0.5, 0.6, True)]
    return cases(signal, (~valid).astype(np.int8), edges, specs)


def nebula_dpdf(rng):
    """NEBULA3D 3D-DeltaPDF: root data (x, y, z), FFT-grid centers in Angstrom, hexagonal cell in lat_* attributes."""
    shape = (8, 8, 6)
    data = rng.normal(0, 1, shape)
    axes = dict(x_axis=(np.arange(8) - 4) * 0.5, y_axis=(np.arange(8) - 4) * 0.5, z_axis=(np.arange(6) - 3) * 0.6)
    with h5py.File(OUT / 'nebula3d_dpdf_small.h5', 'w') as f:
        f.create_dataset('data', data=data, compression='gzip', compression_opts=4)
        for name, values in axes.items():
            f.create_dataset(name, data=values)
        f.attrs.update(q_max=5.0, apodization='gaussian', source_file='synthetic_backfilled.h5')
        for k, v in dict(a=4.0, b=4.0, c=6.0, alpha=90.0, beta=90.0, gamma=120.0).items():
            f.attrs[f'lat_{k}'] = v
    edges = [centers_to_edges(axes[k]) for k in ('z_axis', 'y_axis', 'x_axis')]
    specs = [(2, 0.0, 0.5, False), (0, 0.6, 0.6, False)]   # inversion by index reversal needs a grid symmetric about 0
    return cases(data, np.zeros(shape, np.int8), edges, specs)


def mdhisto_dpdf(rng):
    """3D-DeltaPDF in the Mantid layout, as NEBULA3D writes it: axes x, y, z in Angstrom
    (D2, D1, D0), errors_squared all zero, the cell in the oriented lattice's UB."""
    shape = (8, 8, 6)
    data = rng.normal(0, 1, shape)
    axes = dict(x=(np.arange(8) - 4) * 0.5, y=(np.arange(8) - 4) * 0.5, z=(np.arange(6) - 3) * 0.6)
    params = dict(a=4.0, b=4.0, c=6.0, alpha=90.0, beta=90.0, gamma=120.0)
    with h5py.File(OUT / 'mdhisto_dpdf_small.nxs', 'w') as f:
        entry = f.create_group('MDHistoWorkspace')
        entry.attrs['NX_class'] = 'NXentry'
        group = entry.create_group('data')
        group.attrs['NX_class'] = 'NXdata'
        ds = group.create_dataset('signal', data=data)
        ds.attrs['signal'] = 1
        ds.attrs['axes'] = np.bytes_('D2:D1:D0')
        group.create_dataset('errors_squared', data=np.zeros(shape))
        group.create_dataset('mask', data=np.zeros(shape, np.int8))
        for i, name in enumerate(['z', 'y', 'x']):
            d = group.create_dataset(f'D{i}', data=centers_to_edges(axes[name]))
            d.attrs.update(frame=np.bytes_('General Frame'), long_name=np.bytes_(name), units=np.bytes_('Angstrom'))
        cell = entry.create_group('experiment0/sample/oriented_lattice')
        cell.attrs['NX_class'] = 'NXcrystal'
        cell.create_dataset('orientation_matrix', data=orientation_matrix(**params))
    edges = [centers_to_edges(axes[k]) for k in ('z', 'y', 'x')]
    return cases(data, np.zeros(shape, np.int8), edges, [(2, 0.0, 0.5, False), (1, 0.4, 0.4, False)])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260928)
    expected = dict(mdhisto_small=mdhisto(rng), nxdata_small=nxdata(rng), nebula3d_small=nebula(rng),
                    nebula3d_dpdf_small=nebula_dpdf(rng), mdhisto_dpdf_small=mdhisto_dpdf(rng))
    (OUT / 'expected.json').write_text(json.dumps(expected))
    print('wrote', *sorted(p.name for p in OUT.iterdir()))


if __name__ == '__main__':
    main()
