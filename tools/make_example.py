"""Write examples/demo_hexagonal.nxs, a small synthetic dataset for the demo.

A hexagonal crystal (a = 4.2, c = 6.8 Angstrom) in the Mantid MDHistoWorkspace
layout: Bragg peaks with a two-site structure factor, diffuse rods along L,
partial coverage in wedges (so symmetry averaging has gaps to fill) and
spuriously bright voxels along the coverage edges (for the mask to remove).
Requires numpy and h5py.
"""
from pathlib import Path

import h5py
import numpy as np

OUT = Path(__file__).resolve().parents[1] / 'examples' / 'demo_hexagonal.nxs'
N, HALF = 101, 5.05                      # 101 bins of 0.1 r.l.u. per axis
A, C = 4.2, 6.8


def reciprocal_basis():
    """B with Q = B @ hkl (Cartesian, 1/Angstrom, no 2*pi) and B^T B the reciprocal metric."""
    g = np.array([[A * A, -A * A / 2, 0], [-A * A / 2, A * A, 0], [0, 0, C * C]])
    return np.linalg.cholesky(np.linalg.inv(g)).T


def main():
    rng = np.random.default_rng(7)
    edges = np.linspace(-HALF, HALF, N + 1)
    c = (edges[1:] + edges[:-1]) / 2
    L, K, H = np.meshgrid(c, c, c, indexing='ij')          # storage order (L, K, H)
    B = reciprocal_basis()
    Q = np.stack([H, K, L], -1) @ B.T                      # Cartesian Q = B @ hkl
    q2 = (Q ** 2).sum(-1)

    # Bragg peaks at integer hkl: two sites at (1/3, 2/3, 1/4) and (2/3, 1/3, 3/4).
    signal = np.zeros(H.shape)
    hkl = np.stack(np.meshgrid(*[np.arange(-5, 6)] * 3, indexing='ij'), -1).reshape(-1, 3)
    for h, k, l in hkl:
        f = np.exp(2j * np.pi * (h / 3 + 2 * k / 3 + l / 4)) + np.exp(2j * np.pi * (2 * h / 3 + k / 3 + 3 * l / 4))
        i = abs(f) ** 2 * np.exp(-0.6 * (np.array([h, k, l]) @ B.T @ B @ np.array([h, k, l])))
        if i < 1e-3:
            continue
        d2 = (H - h) ** 2 + (K - k) ** 2 + (L - l) ** 2
        near = d2 < 0.25 ** 2
        signal[near] += 2000 * i * np.exp(-d2[near] / (2 * 0.06 ** 2))
    # Diffuse rods along L through the K points (1/3, 1/3) and equivalents.
    for kh, kk in [(1 / 3, 1 / 3), (2 / 3, -1 / 3), (-1 / 3, 2 / 3), (-1 / 3, -1 / 3), (-2 / 3, 1 / 3), (1 / 3, -2 / 3)]:
        for dh in range(-5, 6):
            for dk in range(-5, 6):
                d2 = (H - kh - dh) ** 2 + (K - kk - dk) ** 2 - (H - kh - dh) * (K - kk - dk)
                signal += 12 * np.exp(-d2 / (2 * 0.05 ** 2)) * (1 + np.cos(np.pi * L) ** 2) / 2
    signal *= np.exp(-0.25 * q2)
    signal += 3 * np.exp(-0.8 * q2)                        # smooth background
    signal = rng.normal(signal, np.sqrt(signal + 0.5) * 0.35)

    # Coverage: five wedges about c*, |Q| < 1.35, and a small beamstop hole.
    phi = np.degrees(np.arctan2(Q[..., 1], Q[..., 0])) % 360
    wedges = [(10, 75), (95, 140), (165, 230), (255, 300), (320, 350)]
    covered = np.zeros(H.shape, bool)
    for lo, hi in wedges:
        covered |= (phi >= lo) & (phi < hi)
    covered &= (q2 < 1.35 ** 2) & (q2 > 0.12 ** 2)
    # Detector-edge artifacts: bright voxels on the coverage boundary.
    inside = covered.copy()
    for axis in range(3):
        for shift in (-1, 1):
            inside &= np.roll(covered, shift, axis)
    edge = covered & ~inside
    signal[edge] *= rng.uniform(4, 12, edge.sum())
    signal[edge] += rng.uniform(40, 200, edge.sum())
    signal[~covered] = np.nan
    signal = np.round(signal, 1).astype(np.float32)        # rounding helps compression

    OUT.parent.mkdir(exist_ok=True)
    with h5py.File(OUT, 'w') as f:
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
    print(f'wrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB), {covered.mean():.1%} covered, {edge.sum()} edge voxels')


if __name__ == '__main__':
    main()
