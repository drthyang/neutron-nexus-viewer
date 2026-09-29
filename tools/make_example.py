"""Write the synthetic example datasets in examples/.

A hexagonal crystal (a = 4.2, c = 6.8 Angstrom) measured at two temperatures,
in the Mantid MDHistoWorkspace layout, telling an order-disorder story:

demo_300K.nxs: Bragg peaks from a two-site structure factor with thermal
diffuse halos, and short-range-order diffuse scattering at the M points:
rods along L (the correlations are mostly in-plane), strongest at even L.

demo_10K.nxs: below the ordering transition the M-point diffuse scattering
condenses into sharp superlattice peaks at even L, the halos fade, and the
smaller Debye-Waller factor strengthens high-Q Bragg peaks.

Both have partial coverage in wedges (so symmetry averaging has gaps to
fill) and spuriously bright voxels along the coverage edges (for the mask to
remove). Requires numpy and h5py.
"""
from pathlib import Path

import h5py
import numpy as np

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'
N, HALF = 101, 5.05                      # 101 bins of 0.1 r.l.u. per axis
A, C = 4.2, 6.8
# The three M points of the hexagonal zone, (1/2, 0), (0, 1/2) and (1/2, -1/2).
M_POINTS = [(0.5, 0.0), (0.0, 0.5), (0.5, -0.5)]


def reciprocal_basis():
    """B with Q = B @ hkl (Cartesian, 1/Angstrom, no 2*pi) and B^T B the reciprocal metric."""
    g = np.array([[A * A, -A * A / 2, 0], [-A * A / 2, A * A, 0], [0, 0, C * C]])
    return np.linalg.cholesky(np.linalg.inv(g)).T


def hex_d2(dh, dk):
    """In-plane squared distance in units of a* (hexagonal metric, gamma* = 60 degrees)."""
    return dh * dh + dk * dk + dh * dk


def near(h, k, l, r_hk, r_l):
    """Storage-order (L, K, H) slices of the bins within r_hk, r_l r.l.u. of (h, k, l)."""
    def axis(x, r):
        return slice(max(0, int(np.floor((x - r + HALF) / 0.1))), min(N, int(np.ceil((x + r + HALF) / 0.1)) + 1))
    return axis(l, r_l), axis(k, r_hk), axis(h, r_hk)


def write(out, *, temperature, ordered, rotation):
    rng = np.random.default_rng(7 if ordered else 11)
    edges = np.linspace(-HALF, HALF, N + 1)
    c = (edges[1:] + edges[:-1]) / 2
    L, K, H = np.meshgrid(c, c, c, indexing='ij')          # storage order (L, K, H)
    B = reciprocal_basis()
    Q = np.stack([H, K, L], -1) @ B.T                      # Cartesian Q = B @ hkl
    q2 = (Q ** 2).sum(-1)
    debye_waller = 0.12 if ordered else 0.4                 # exp(-u q^2), larger when hot
    signal = np.zeros(H.shape)

    # Bragg peaks at integer hkl, two sites at (1/3, 2/3, 1/4) and (2/3, 1/3, 3/4),
    # each with a thermal diffuse halo ~ T q^2 / (dq^2 + kappa^2).
    hkl = np.stack(np.meshgrid(*[np.arange(-5, 6)] * 3, indexing='ij'), -1).reshape(-1, 3)
    halo = 0.02 * temperature / 300
    for h, k, l in hkl:
        f = np.exp(2j * np.pi * (h / 3 + 2 * k / 3 + l / 4)) + np.exp(2j * np.pi * (2 * h / 3 + k / 3 + 3 * l / 4))
        qg2 = np.array([h, k, l]) @ B.T @ B @ np.array([h, k, l])
        i = abs(f) ** 2 * np.exp(-(0.35 + debye_waller) * qg2)
        if i < 1e-3:
            continue
        box = near(h, k, l, 0.8, 1.4)
        d2 = hex_d2(H[box] - h, K[box] - k) * B[0, 0] ** 2 + ((L[box] - l) / C) ** 2   # 1/A^2
        signal[box] += i * (1500 * np.exp(-d2 / (2 * 0.016 ** 2))
                            + 3 * halo * qg2 / (d2 + 0.012 ** 2) * np.exp(-d2 / (2 * 0.07 ** 2)))

    # M-point order: short-range (hot) as broad in-plane Lorentzian-squared rods along
    # L, long-range (cold) as sharp peaks; both strongest at even L.
    even_l = np.cos(np.pi * L / 2) ** 2
    rods = np.zeros(H.shape[1:])                           # (K, H) plane
    width = 0.07 if ordered else 0.15
    for dh in range(-6, 6):
        for dk in range(-6, 6):
            for mh, mk in M_POINTS:
                d2 = hex_d2(H[0] - dh - mh, K[0] - dk - mk)
                rods += 1 / (1 + d2 / width ** 2) ** 2
    signal += (4 if ordered else 60) * rods[None] * (0.15 + 0.85 * even_l) * np.exp(-0.5 * q2)
    if ordered:
        for dh in range(-6, 6):
            for dk in range(-6, 6):
                for mh, mk in M_POINTS:
                    ch, ck = dh + mh, dk + mk
                    if abs(ch) > HALF or abs(ck) > HALF:
                        continue
                    for l in range(-4, 5, 2):
                        box = near(ch, ck, l, 0.4, 0.6)
                        d2 = hex_d2(H[box] - ch, K[box] - ck) * B[0, 0] ** 2 + ((L[box] - l) / C) ** 2
                        signal[box] += 300 * np.exp(-d2 / (2 * 0.016 ** 2)) * np.exp(-0.5 * q2[box])

    signal += 2 + 5 * np.exp(-0.8 * q2)                    # smooth background
    signal = rng.normal(signal, np.sqrt(signal + 0.5) * 0.3)

    # Coverage: four wedges about c* (rotated per dataset), 0.12 < |Q| < 1.35.
    phi = (np.degrees(np.arctan2(Q[..., 1], Q[..., 0])) - rotation) % 360
    covered = np.zeros(H.shape, bool)
    for lo, hi in [(5, 80), (95, 170), (185, 260), (275, 350)]:
        covered |= (phi >= lo) & (phi < hi)
    covered &= (q2 < 1.35 ** 2) & (q2 > 0.12 ** 2)
    # Detector-edge artifacts: bright voxels on the coverage boundary.
    inside = covered.copy()
    for axis in range(3):
        for shift in (-1, 1):
            inside &= np.roll(covered, shift, axis)
    edge = covered & ~inside
    signal[edge] *= rng.uniform(3, 8, edge.sum())
    signal[edge] += rng.uniform(30, 120, edge.sum())
    signal[~covered] = np.nan
    signal = np.round(signal, 1).astype(np.float32)        # rounding helps compression

    out.parent.mkdir(exist_ok=True)
    with h5py.File(out, 'w') as f:
        entry = f.create_group('MDHistoWorkspace')
        entry.attrs['NX_class'] = 'NXentry'
        data = entry.create_group('data')
        data.attrs['NX_class'] = 'NXdata'
        kw = dict(chunks=(1, N, N), compression='gzip', compression_opts=9, shuffle=True)
        ds = data.create_dataset('signal', data=signal, **kw)
        ds.attrs['signal'] = 1
        ds.attrs['axes'] = np.bytes_('D2:D1:D0')
        data.create_dataset('mask', data=np.zeros(signal.shape, np.int8), **kw)
        for i, name in enumerate(['[H,0,0]', '[0,K,0]', '[0,0,L]']):
            d = data.create_dataset(f'D{i}', data=edges)
            d.attrs.update(frame=np.bytes_('HKL'), long_name=np.bytes_(name), units=np.bytes_('r.l.u.'))
        lattice = entry.create_group('experiment0/sample/oriented_lattice')
        lattice.attrs['NX_class'] = 'NXcrystal'
        for k, v in dict(a=A, b=A, c=C, alpha=90.0, beta=90.0, gamma=120.0).items():
            lattice.create_dataset(f'unit_cell_{k}', data=[v])
        lattice.create_dataset('orientation_matrix', data=B)
    print(f'wrote {out} ({out.stat().st_size / 1e6:.1f} MB), {covered.mean():.1%} covered, {edge.sum()} edge voxels')


def main():
    write(EXAMPLES / 'demo_300K.nxs', temperature=300, ordered=False, rotation=0)
    write(EXAMPLES / 'demo_10K.nxs', temperature=10, ordered=True, rotation=8)


if __name__ == '__main__':
    main()
