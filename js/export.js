// Export of the processed volume as an input file for NEBULA3D's 3D-ΔPDF pipeline
// (https://github.com/drthyang/nebula3d).
//
// NEBULA3D reads an HDF5 group /entry with `data` in (H, K, L) C order (H slowest),
// `mask` (1 = valid), `h_axis`, `k_axis`, `l_axis` at bin centres and `ub_matrix`
// with the 2π factor (|Q| = |UB·hkl|). It does not symmetrize, and its ΔPDF puts
// Q = 0 at index n//2 of every axis, so the grid is padded to be symmetric about
// 0. Every voxel of the padded grid gets the equal-weight mean of the measured,
// unmasked voxels among its symmetry equivalents, so symmetry also fills the
// padding; voxels with none are written as 0 with mask 0, which NEBULA3D backfills.
//
// Conventions as in slab.js: `shape` is the storage (C-order) shape, display
// dimension d lives on storage axis 2 - d, and symmetry maps are integer affine
// maps i' = M i + t on display-dimension bin indices (from indexMaps()).

import { reciprocalMetric } from './nexus.js';

const near = (x, tol = 1e-3) => Math.abs(x - Math.round(x)) < tol;

/**
 * The exported grid: which display dimension holds H, K and L, and how each is
 * padded to be symmetric about 0. Throws an Error that says why the data cannot
 * be exported (no unit cell, projected or non-HKL axes, non-uniform bins, or a
 * grid with neither a bin centre nor a bin edge at 0).
 */
export function exportPlan(dims, lattice) {
  if (!lattice) throw new Error('Export needs the unit cell (a UB matrix or unit_cell_* in the file): NEBULA3D computes |Q| from it.');
  const order = [0, 1, 2].map((a) => {
    const d = dims.findIndex((dim) => dim.basis && dim.basis.vec.every((v, j) => v === (j === a ? 1 : 0)));
    if (d < 0) throw new Error(`Export needs plain H, K and L axes; this file has ${dims.map((x) => x.label).join(', ')}.`);
    return d;
  });
  const lo = [0, 0, 0], size = [0, 0, 0], source = [0, 0, 0], paddedDims = dims.slice();
  dims.forEach((dim, d) => {
    const e = dim.edges, n = e.length - 1, w = (e[n] - e[0]) / n;
    for (let k = 0; k < n; k++) {
      if (Math.abs(e[k + 1] - e[k] - w) > 1e-4 * w) throw new Error(`Axis ${dim.label} has non-uniform bins; export needs uniform bins.`);
    }
    const centre = -e[0] / w - 0.5, edge = -e[0] / w;
    if (near(centre)) {
      // A bin centred at 0: an odd grid with that bin in the middle.
      const i0 = Math.round(centre), half = Math.max(i0, n - 1 - i0);
      lo[d] = i0 - half;
      size[d] = 2 * half + 1;
    } else if (near(edge)) {
      // A bin edge at 0: an even grid split there.
      const j0 = Math.round(edge), half = Math.max(j0, n - j0);
      lo[d] = j0 - half;
      size[d] = 2 * half;
    } else {
      throw new Error(`Axis ${dim.label} has neither a bin centre nor a bin edge at 0; NEBULA3D needs a grid symmetric about Q = 0.`);
    }
    if (size[d] > 3 * n) throw new Error(`Axis ${dim.label} (${e[0].toFixed(2)} to ${e[n].toFixed(2)}) is too far from symmetric about 0 to pad.`);
    source[d] = n;
    paddedDims[d] = { ...dim, edges: Float64Array.from({ length: size[d] + 1 }, (_, k) => e[0] + (lo[d] + k) * w) };
  });
  const centers = order.map((d) => {
    const e = paddedDims[d].edges;
    return Float64Array.from({ length: size[d] }, (_, k) => (e[k] + e[k + 1]) / 2);
  });
  return {
    order, lo, size, source, paddedDims, centers,
    shape: order.map((d) => size[d]),
    padded: size.some((s, d) => s !== source[d]),
    ub: nebulaUB(lattice),
  };
}

/**
 * NEBULA3D's UB (row-major, with 2π): the file's orientation matrix (Mantid,
 * without 2π) times 2π, or without one, the Cholesky factor of the reciprocal
 * metric. Only the metric (UB)ᵀUB matters for |Q| and the real-space axes.
 */
export function nebulaUB(lattice) {
  let ub = lattice.ub;
  if (!ub) {
    const G = reciprocalMetric(lattice);
    const b00 = Math.sqrt(G[0][0]), b01 = G[0][1] / b00, b02 = G[0][2] / b00;
    const b11 = Math.sqrt(G[1][1] - b01 * b01), b12 = (G[1][2] - b01 * b02) / b11;
    const b22 = Math.sqrt(G[2][2] - b02 * b02 - b12 * b12);
    ub = [b00, b01, b02, 0, b11, b12, 0, 0, b22];
  }
  return ub.map((x) => 2 * Math.PI * x);
}

/**
 * The symmetrized volume on the padded grid of `plan`, in (H, K, L) C order.
 * `maps` are the group's index maps on plan.paddedDims (identity included); a
 * user `mask` (nonzero = removed, storage order) excludes voxels. Each orbit is
 * visited once, from its smallest padded flat index; every member gets the mean
 * of the orbit's valid source voxels (repeated images scale sum and count alike).
 * Returns { data: Float32Array, valid: Uint8Array, stats }.
 */
export function symmetrizeForExport(volume, shape, plan, maps, mask = null, onProgress = null) {
  const n = [shape[2], shape[1], shape[0]];
  const [P0, P1, P2] = plan.size, [lo0, lo1, lo2] = plan.lo;
  // Output strides per display dimension: C order over (H, K, L).
  const [nH, nK, nL] = plan.shape, outStride = [nK * nL, nL, 1];
  const os = [0, 0, 0];
  plan.order.forEach((d, a) => { os[d] = outStride[a]; });
  const total = P0 * P1 * P2, data = new Float32Array(total), valid = new Uint8Array(total);
  const G = maps.length, M = new Int32Array(9 * G), T = new Int32Array(3 * G);
  maps.forEach((m, g) => { M.set(m.M, 9 * g); T.set(m.t, 3 * g); });
  const src = new Int32Array(G), out = new Int32Array(G);
  let validCount = 0, measured = 0;

  for (let p2 = 0; p2 < P2; p2++) {
    for (let p1 = 0; p1 < P1; p1++) {
      for (let p0 = 0; p0 < P0; p0++) {
        const flat = p0 + P0 * (p1 + P1 * p2);
        let k = 0, rep = true;
        for (let g = 0; g < G; g++) {
          const o = 9 * g, t = 3 * g;
          const q0 = M[o] * p0 + M[o + 1] * p1 + M[o + 2] * p2 + T[t];
          if (q0 < 0 || q0 >= P0) continue;
          const q1 = M[o + 3] * p0 + M[o + 4] * p1 + M[o + 5] * p2 + T[t + 1];
          if (q1 < 0 || q1 >= P1) continue;
          const q2 = M[o + 6] * p0 + M[o + 7] * p1 + M[o + 8] * p2 + T[t + 2];
          if (q2 < 0 || q2 >= P2) continue;
          if (q0 + P0 * (q1 + P1 * q2) < flat) { rep = false; break; }
          const s0 = q0 + lo0, s1 = q1 + lo1, s2 = q2 + lo2;
          src[k] = s0 >= 0 && s0 < n[0] && s1 >= 0 && s1 < n[1] && s2 >= 0 && s2 < n[2] ? s0 + n[0] * (s1 + n[1] * s2) : -1;
          out[k] = q0 * os[0] + q1 * os[1] + q2 * os[2];
          k++;
        }
        if (!rep) continue;
        let sum = 0, count = 0;
        for (let i = 0; i < k; i++) {
          const s = src[i];
          if (s < 0) continue;
          const v = volume[s];
          if (v === v && !(mask && mask[s])) { sum += v; count++; }
        }
        if (!count) continue;
        const mean = sum / count;
        for (let i = 0; i < k; i++) {
          if (!valid[out[i]]) validCount++;
          data[out[i]] = mean;
          valid[out[i]] = 1;
        }
      }
    }
    onProgress?.((p2 + 1) / P2);
  }
  for (let i = 0; i < volume.length; i++) if (volume[i] === volume[i] && !(mask && mask[i])) measured++;
  return { data, valid, stats: { total, valid: validCount, measured, shape: [nH, nK, nL] } };
}

/**
 * Write a NEBULA3D input file at `path` with h5wasm: /entry/{data, mask,
 * h_axis, k_axis, l_axis, ub_matrix}, gzip-compressed, with provenance `attrs`
 * (strings) on /entry.
 */
export function writeNebulaFile(h5wasm, path, plan, result, attrs = {}) {
  const f = new h5wasm.File(path, 'w');
  try {
    const entry = f.create_group('entry');
    const shape = plan.shape, chunks = [1, shape[1], shape[2]];
    entry.create_dataset({ name: 'data', data: result.data, shape, chunks, compression: 'gzip', compression_opts: 1 });
    entry.create_dataset({ name: 'mask', data: result.valid, shape, chunks, compression: 'gzip', compression_opts: 1 });
    ['h_axis', 'k_axis', 'l_axis'].forEach((name, a) => entry.create_dataset({ name, data: plan.centers[a] }));
    entry.create_dataset({ name: 'ub_matrix', data: Float64Array.from(plan.ub), shape: [3, 3] });
    for (const [key, value] of Object.entries({ instrument: '', ...attrs })) entry.create_attribute(key, String(value));
  } finally {
    f.close();
  }
}
