// Reindex an HKL volume with another UB matrix (pure functions): read an ISAW UB file,
// relate the two indexings, plan the new grid and resample the volume onto it.
//
// Both UBs map indices to the same Q in the sample frame (Mantid's convention, no 2π):
// Q = UB_old·h = UB_new·h′, so h′ = T·h with T = UB_new⁻¹·UB_old. A display coordinate u
// of the volume (one per axis, along the axis's basis vector, so h = W·u) is at
// h′ = T·W·u, and a new voxel centre h′ reads the old volume at u = (T·W)⁻¹·h′.

import { isHKL } from './nexus.js';
import { det3, inv3, mul, transpose } from './rigaku-geometry.js';

/**
 * An ISAW UB file (*.mat, as Mantid's SaveIsawUB writes it): three rows of the transposed
 * UB in ISAW's frame (x along the beam, z up), then a b c α β γ and the cell volume.
 * Returns the UB in Mantid's frame (z along the beam, y up), row-major and without 2π,
 * as Mantid's LoadIsawUB reads it, and the cell the file lists (or null).
 */
export function parseIsawUB(text) {
  const rows = [];
  for (const line of String(text).split(/\r?\n|\r/)) {
    const words = line.trim().split(/\s+/).filter(Boolean);
    if (!words.length) continue;
    const v = words.map(Number);
    if (!v.every(Number.isFinite)) break; // the comment lines after the numbers
    rows.push(v);
  }
  if (rows.length < 3 || rows.slice(0, 3).some((r) => r.length < 3)) {
    throw new Error('This is not an ISAW UB file: it should start with three rows of three numbers (the transposed UB).');
  }
  const F = rows.slice(0, 3).map((r) => r.slice(0, 3));
  // File row j = column j of the ISAW UB; Mantid row i = ISAW row (i + 1) mod 3.
  const ub = [];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) ub.push(F[j][(i + 1) % 3]);
  if (!(Math.abs(det3(ub)) > 1e-12)) throw new Error('The UB matrix in this file is singular.');
  const listed = rows.slice(3).find((r) => r.length >= 6 && r.slice(0, 3).every((x) => x > 0) && r.slice(3, 6).every((x) => x > 0 && x < 180));
  const cell = listed ? { a: listed[0], b: listed[1], c: listed[2], alpha: listed[3], beta: listed[4], gamma: listed[5] } : null;
  return { ub, cell };
}

/** The display→hkl matrix W: column d is display axis d's basis vector (h = W·u). */
export function basisMatrix(dims) {
  const W = new Array(9);
  for (let i = 0; i < 3; i++) for (let d = 0; d < 3; d++) W[3 * i + d] = dims[d].basis.vec[i];
  return W;
}

/** Each display axis as { e0, step, n } (first edge, bin width, bins); throws unless the bins are uniform. */
export function uniformAxes(dims) {
  return dims.map((d) => {
    const e = d.edges, n = e.length - 1, step = (e[n] - e[0]) / n;
    for (let k = 1; k <= n; k++) {
      if (Math.abs(e[k] - e[0] - k * step) > 1e-4 * step) throw new Error(`axis ${d.label} has bins of unequal width; reindexing needs uniform bins.`);
    }
    return { e0: e[0], step, n };
  });
}

/** Why a dataset cannot be reindexed (a sentence), or '' when it can. */
export function reindexProblem(info) {
  if (!info) return 'Open a dataset first.';
  if (!info.dims.every((d) => d.basis && isHKL(d))) return 'Reindexing needs axes in r.l.u. along HKL directions.';
  if (!info.lattice?.ub) return 'This file has no UB matrix (orientation_matrix), so the orientation its HKL axes were made with is unknown.';
  try {
    uniformAxes(info.dims);
  } catch (err) {
    return `Reindexing needs uniform bins: ${err.message}`;
  }
  return '';
}

/** h′ = T·h for indices h of UB `ubOld` and h′ of UB `ubNew` (both mapping to the same Q). */
export const reindexMatrix = (ubOld, ubNew) => mul(inv3(ubNew), ubOld);

/**
 * How T relates the two indexings: the nearest integer matrix N and the largest
 * deviation from it, and, when N is a change of setting (nonsingular, right-handed), the
 * rotation angle (°) between the orientations once the old UB is put in the new
 * setting (UB_old·N⁻¹), from the polar decomposition of UB_new·N·UB_old⁻¹.
 */
export function describeTransform(T, ubOld, ubNew) {
  const N = T.map((x) => Math.round(x) || 0);
  const deviation = Math.max(...T.map((x, i) => Math.abs(x - N[i])));
  let angle = null;
  if (det3(N) > 0.5) {
    let R = mul(mul(ubNew, N), inv3(ubOld));
    if (det3(R) > 0) {
      for (let k = 0; k < 30; k++) {
        const Rt = transpose(inv3(R));
        R = R.map((v, i) => 0.5 * (v + Rt[i]));
      }
      angle = Math.acos(Math.max(-1, Math.min(1, (R[0] + R[4] + R[8] - 1) / 2))) * 180 / Math.PI;
    }
  }
  return { N, deviation, angle, det: det3(T) };
}

/** The matrices of a reindexing: W (h = W·u), TW (h′ = TW·u) and M = (TW)⁻¹ (u = M·h′), and the uniform axes. */
export function reindexGeometry(dims, T) {
  const W = basisMatrix(dims), TW = mul(T, W);
  return { W, TW, M: inv3(TW), axes: uniformAxes(dims) };
}

/**
 * The bounding box, in new indices, of the measured voxels (finite values of `volume`,
 * in the order of `shape`, display axis d on volume axis 2 - d), padded by half a voxel's
 * extent along each new axis: { min, max } or null when nothing is measured.
 */
export function measuredBox(volume, shape, axes, TW) {
  const [n0, n1, n2] = shape;
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  const u = axes.map((a) => Float64Array.from({ length: a.n }, (_, k) => a.e0 + (k + 0.5) * a.step));
  for (let i0 = 0, v = 0; i0 < n0; i0++) {
    const u2 = u[2][i0];
    for (let i1 = 0; i1 < n1; i1++) {
      const u1 = u[1][i1];
      for (let i2 = 0; i2 < n2; i2++, v++) {
        if (volume[v] !== volume[v]) continue;
        const u0 = u[0][i2];
        for (let i = 0; i < 3; i++) {
          const h = TW[3 * i] * u0 + TW[3 * i + 1] * u1 + TW[3 * i + 2] * u2;
          if (h < min[i]) min[i] = h;
          if (h > max[i]) max[i] = h;
        }
      }
    }
  }
  if (!(min[0] <= max[0])) return null;
  const pad = [0, 1, 2].map((i) => 0.5 * axes.reduce((s, a, d) => s + Math.abs(TW[3 * i + d]) * a.step, 0));
  return { min: min.map((x, i) => x - pad[i]), max: max.map((x, i) => x + pad[i]) };
}

// Voxel sizes 1/n, so that integer indices fall on voxel centres.
const NICE_N = [1, 2, 4, 5, 8, 10, 16, 20, 25, 40, 50, 80, 100, 125, 200, 250, 400, 500, 1000];

/**
 * A voxel size per new axis that keeps the old grid's resolution: a step of s along new
 * axis i moves s·M[:, i] in display coordinates, s/‖(M[d][i]/Δ_d)_d‖ is one old voxel,
 * rounded to the nearest 1/n of NICE_N or an old bin width that divides 1.
 */
export function defaultSteps(M, axes) {
  const candidates = NICE_N.map((n) => 1 / n);
  for (const a of axes) {
    for (const s of [a.step, 2 * a.step, a.step / 2]) if (Math.abs(1 / s - Math.round(1 / s)) < 1e-6 * (1 / s)) candidates.push(1 / Math.round(1 / s));
  }
  return [0, 1, 2].map((i) => {
    const raw = 1 / Math.hypot(...axes.map((a, d) => M[3 * d + i] / a.step));
    let best = candidates[0];
    for (const c of candidates) {
      const better = Math.abs(Math.log(c / raw)) - Math.abs(Math.log(best / raw));
      if (better < -1e-9 || (Math.abs(better) <= 1e-9 && c < best)) best = c;
    }
    return best;
  });
}

/** Sub-samples per axis: about two per old voxel the new voxel spans (1 to 5). */
export function defaultSubsamples(M, axes, steps) {
  const extent = Math.max(...axes.map((a, d) => steps.reduce((s, si, i) => s + Math.abs(M[3 * d + i]) * si, 0) / a.step));
  return Math.min(5, Math.max(1, Math.round(2 * extent)));
}

/**
 * The new grid over `box` with voxel sizes `steps`: voxel centres at integer multiples of
 * the step, so integer indices are centres when 1/step is an integer. Returns
 * { min (first centre per axis H, K, L), steps, shape [nH, nK, nL], voxels, edges }.
 */
export function reindexGrid(box, steps) {
  // A box edge within 0.001 voxel of a boundary (a UB read with 8 decimals moves it by
  // about 1e-6) does not open a voxel beyond it.
  const k0 = box.min.map((x, i) => Math.floor(x / steps[i] + 0.5 + 1e-3));
  const k1 = box.max.map((x, i) => Math.ceil(x / steps[i] - 0.5 - 1e-3));
  const shape = k0.map((k, i) => k1[i] - k + 1);
  const edges = shape.map((n, i) => Array.from({ length: n + 1 }, (_, k) => Number(((k0[i] + k - 0.5) * steps[i]).toPrecision(12))));
  return { min: k0.map((k, i) => k * steps[i]), steps: [...steps], shape, voxels: shape[0] * shape[1] * shape[2], edges };
}

/**
 * Resample `volume` (float32 in the order of `shape`, NaN where unmeasured) onto `grid`:
 * each new voxel is the equal-weight mean of nsub³ sub-samples spread over it, each
 * taking the old voxel it falls in (unmeasured ones are left out). With `variance`, σ² of
 * the mean is Σ_k n_k² σ_k² / N², where n_k of the N measured sub-samples fall in old
 * voxel k (a voxel sampled twice counts as one correlated value). Returns, in Mantid's
 * storage order (H fastest, L slowest), signal (NaN where no sub-sample was measured),
 * errors2 (0 without variance) and events (measured sub-samples), and the covered voxels.
 */
export function reindexVolume(volume, variance, shape, axes, M, grid, nsub, onProgress = () => {}) {
  const [nH, nK, nL] = grid.shape, n = nH * nK * nL;
  const [n0, n1, n2] = shape; // n2 bins along display axis 0, n0 along display axis 2
  // Fractional old voxel indices along display axis d: f_d = A[d]·h′ + o_d.
  const A = M.map((x, k) => x / axes[Math.floor(k / 3)].step), o = axes.map((a) => -a.e0 / a.step);
  const P = nsub ** 3, S = new Float64Array(3 * P), s = grid.steps;
  for (let p = 0; p < P; p++) {
    const q = [p % nsub, Math.floor(p / nsub) % nsub, Math.floor(p / nsub / nsub)];
    const delta = q.map((x, i) => ((x + 0.5) / nsub - 0.5) * s[i]);
    for (let d = 0; d < 3; d++) S[3 * p + d] = A[3 * d] * delta[0] + A[3 * d + 1] * delta[1] + A[3 * d + 2] * delta[2];
  }
  let signal, errors2, events;
  try {
    signal = new Float64Array(n);
    errors2 = new Float64Array(n);
    events = new Float64Array(n);
  } catch {
    throw new Error(`Could not allocate ${(n * 24 / 1e9).toFixed(2)} GB for ${n.toLocaleString()} voxels; use a larger voxel.`);
  }
  const ids = new Int32Array(P), counts = new Int32Array(P);
  const [a0, a1, a2] = [A[0] * s[0], A[3] * s[0], A[6] * s[0]];
  let covered = 0;
  for (let c = 0, out = 0; c < nL; c++) {
    const hL = grid.min[2] + c * s[2];
    for (let b = 0; b < nK; b++) {
      const hK = grid.min[1] + b * s[1], hH = grid.min[0];
      const f0 = A[0] * hH + A[1] * hK + A[2] * hL + o[0];
      const f1 = A[3] * hH + A[4] * hK + A[5] * hL + o[1];
      const f2 = A[6] * hH + A[7] * hK + A[8] * hL + o[2];
      for (let a = 0; a < nH; a++, out++) {
        const g0 = f0 + a * a0, g1 = f1 + a * a1, g2 = f2 + a * a2;
        let sum = 0, N = 0, distinct = 0;
        for (let p = 0, t = 0; p < P; p++, t += 3) {
          const j0 = Math.floor(g0 + S[t]), j1 = Math.floor(g1 + S[t + 1]), j2 = Math.floor(g2 + S[t + 2]);
          if (j0 < 0 || j0 >= n2 || j1 < 0 || j1 >= n1 || j2 < 0 || j2 >= n0) continue;
          const idx = j0 + n2 * (j1 + n1 * j2), v = volume[idx];
          if (v !== v) continue;
          sum += v;
          N++;
          if (variance) {
            let k = 0;
            while (k < distinct && ids[k] !== idx) k++;
            if (k === distinct) { ids[distinct] = idx; counts[distinct++] = 0; }
            counts[k]++;
          }
        }
        if (N) {
          signal[out] = sum / N;
          events[out] = N;
          if (variance) {
            let e = 0;
            for (let k = 0; k < distinct; k++) e += counts[k] * counts[k] * variance[ids[k]];
            errors2[out] = e / (N * N);
          }
          covered++;
        } else {
          signal[out] = NaN;
          errors2[out] = NaN;
        }
      }
    }
    onProgress((c + 1) / nL);
  }
  return { signal, errors2, events, covered };
}
