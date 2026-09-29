// Powder average I(|Q|): the 3-D volume reduced to a 1-D curve.
//
// A voxel holds the intensity per unit reciprocal volume in its cell, so the
// average over a shell Q ≤ |Q| < Q + ΔQ is the volume-weighted mean of the
// voxels in it, over the voxels that have data:
//
//   I(Q) = Σ_v f_v I_v / Σ_v f_v,   f_v = the part of voxel v inside the shell.
//
// Unmeasured and masked voxels are left out of both sums rather than counted
// as zero, so incomplete coverage costs precision, not intensity. With a
// symmetry group, an orbit's intensity is the equal-weight mean of its measured
// members (as in the slices), and every distinct member of the orbit counts
// with its own |Q| and volume, measured or not, inside the grid or beyond it.
// Orbits are thus weighted by their multiplicity, not by how many of their
// members happened to be measured, which would bias the average toward the
// directions that were measured best.
//
// Each voxel is split into s³ sub-cells binned by their own |Q|, so a voxel is
// shared between the shells it overlaps (s = 1 bins voxel centres). With
// per-voxel variances (Mantid's errors_squared), the members of an orbit share
// one mean and are fully correlated, so
//
//   σ²(Q) = Σ_o c_o² σ_o² / (Σ_o c_o)²,   σ_o² = Σ_m σ_m² / n_o²,
//
// with c_o the orbit's weight in the shell and the sum over its n_o measured
// members m.
//
// Conventions as in slab.js: `shape` is the storage (C-order) shape, display
// dimension d lives on storage axis 2 - d, and symmetry maps are integer affine
// maps i' = M i + t on display-dimension bin indices (from indexMaps()).

import { isHKL, reciprocalMetric } from './nexus.js';

const MAX_SHELLS = 20000;
const INVERSE_ANGSTROM = /^(?:(?:a|å|ang(?:strom)?s?)\s*(?:\^?\s*-\s*1|⁻¹)|1\s*\/\s*(?:a|å|ang(?:strom)?s?)|inverse\s+ang(?:strom)?s?)$/i;

/**
 * The metric of display coordinates in |Q|: |Q|² = xᵀ G x in Å⁻², with 2π.
 * HKL axes use the reciprocal metric of the cell and the axes' basis vectors;
 * axes in Å⁻¹ (Mantid's Q frames) are taken as Cartesian. Throws an Error that
 * says why |Q| cannot be computed otherwise.
 */
export function qMetric(dims, lattice) {
  if (dims.every((d) => d.basis && isHKL(d))) {
    if (!lattice) throw new Error('I(Q) needs the unit cell (a UB matrix or unit_cell_* in the file) to convert H, K, L to |Q|.');
    const Gs = reciprocalMetric(lattice), k = 4 * Math.PI * Math.PI;
    const dot = (u, v) => u.reduce((s, ui, i) => s + ui * v.reduce((t, vj, j) => t + Gs[i][j] * vj, 0), 0);
    return { G: dims.flatMap((a) => dims.map((b) => k * dot(a.basis.vec, b.basis.vec))), frame: 'HKL' };
  }
  if (dims.every((d) => !isHKL(d) && INVERSE_ANGSTROM.test(d.units.trim()))) return { G: [1, 0, 0, 0, 1, 0, 0, 0, 1], frame: 'Q' };
  throw new Error(`I(Q) needs H, K, L axes and a unit cell, or Q axes in Å⁻¹; this file has ${dims.map((d) => d.units ? `${d.label} (${d.units})` : d.label).join(', ')}.`);
}

/**
 * The shells and sampling for powderAverage(): shells of width `dq` from 0 to
 * `qmax` (Å⁻¹), and voxels split into `split`³ sub-cells. By default ΔQ is the
 * shortest bin step in |Q| and the shells reach the farthest grid corner (`top`).
 */
export function powderPlan(dims, lattice, { dq = null, qmax = null, split = 2 } = {}) {
  const { G, frame } = qMetric(dims, lattice);
  const n = dims.map((d) => d.edges.length - 1);
  const w = dims.map((d, i) => (d.edges[n[i]] - d.edges[0]) / n[i]);
  dims.forEach((d, i) => {
    for (let k = 0; k < n[i]; k++) {
      if (Math.abs(d.edges[k + 1] - d.edges[k] - w[i]) > 1e-4 * w[i]) throw new Error(`Axis ${d.label} has non-uniform bins; I(Q) needs uniform bins.`);
    }
  });
  const qOf = (x) => Math.sqrt(Math.max(0, x.reduce((s, xi, i) => s + xi * x.reduce((t, xj, j) => t + G[3 * i + j] * xj, 0), 0)));
  // |Q| is convex, so its largest value on the grid is at a corner.
  const top = Math.max(...Array.from({ length: 8 }, (_, c) => qOf(dims.map((d, i) => d.edges[(c >> i) & 1 ? n[i] : 0]))));
  const step = Math.min(...w.map((wi, i) => wi * Math.sqrt(G[4 * i])));
  const width = dq ?? Number(step.toPrecision(1));
  const reach = qmax ?? top;
  if (!(width > 0) || !Number.isFinite(width)) throw new Error('ΔQ must be positive.');
  if (!(reach > 0)) throw new Error('Q max must be positive.');
  const count = Math.max(1, Math.ceil(reach / width - 1e-9));
  if (count > MAX_SHELLS) throw new Error(`ΔQ = ${width} Å⁻¹ gives ${count} shells up to ${reach.toFixed(2)} Å⁻¹; use at most ${MAX_SHELLS}.`);
  if (!(Number.isInteger(split) && split >= 1 && split <= 4)) throw new Error('Split voxels into 1 to 4 sub-cells per axis.');
  return {
    G, frame, dq: width, qmax: count * width, count, split, top, step,
    c0: dims.map((d, i) => d.edges[0] + w[i] / 2), w,
    // Volume of one voxel in Å⁻³.
    voxel: Math.sqrt(Math.max(0, det3(G))) * w[0] * w[1] * w[2],
  };
}

/**
 * I(Q) of the volume on the shells of `plan` (see the top of this file).
 * `maps` are the group's index maps on the volume's grid (identity included), a
 * user `mask` (nonzero = removed) and NaN voxels are excluded, and `variance`
 * (σ² per voxel, storage order, or null) gives the uncertainties. Returns per
 * shell the mean `intensity` and its `sigma` (NaN without data or without
 * variances), `voxels` (the number of voxels with data, fractional with split
 * voxels) and `coverage` (the part of the shell's volume that has data).
 */
export function powderAverage(volume, shape, plan, maps, mask = null, variance = null, onProgress = null) {
  const n0 = shape[2], n1 = shape[1], n2 = shape[0], s1 = n0, s2 = n0 * n1;
  const { G, c0, w, dq, count, split } = plan;
  // Sub-cell offsets δ from the voxel centre along each axis (doubled, for
  // the cross term 2δᵀGx), and the constant δᵀGδ of each sub-cell.
  const S = split ** 3, off = (a, d) => ((a + 0.5) / split - 0.5) * w[d];
  const ex = Float64Array.from({ length: split }, (_, a) => 2 * off(a, 0));
  const ey = Float64Array.from({ length: split }, (_, a) => 2 * off(a, 1));
  const ez = Float64Array.from({ length: split }, (_, a) => 2 * off(a, 2));
  const dd = new Float64Array(S);
  for (let a = 0, j = 0; a < split; a++) {
    for (let b = 0; b < split; b++) {
      for (let c = 0; c < split; c++, j++) {
        const d = [off(a, 0), off(b, 1), off(c, 2)];
        dd[j] = d.reduce((s, di, i) => s + di * d.reduce((t, dk, k) => t + G[3 * i + k] * dk, 0), 0);
      }
    }
  }
  const Ng = maps.length, M = new Int32Array(9 * Ng), T = new Int32Array(3 * Ng);
  maps.forEach((m, g) => { M.set(m.M, 9 * g); T.set(m.t, 3 * g); });
  const img = new Int32Array(3 * Ng);
  // Voxels already visited as members of an earlier orbit (one bit each).
  const seen = Ng > 1 ? new Int32Array(Math.ceil(volume.length / 32)) : null;
  const weight = new Float64Array(count), sum = new Float64Array(count), spread = new Float64Array(count);
  // Shell bounds in |Q|², and the shells one orbit falls in with its sub-cell count in each.
  const bound = Float64Array.from({ length: count + 1 }, (_, b) => (b * dq) ** 2);
  const shells = new Int32Array(Ng * S), hits = new Int32Array(Ng * S);
  const inv = 1 / dq;
  let orbits = 0;

  for (let p2 = 0; p2 < n2; p2++) {
    for (let p1 = 0; p1 < n1; p1++) {
      for (let p0 = 0; p0 < n0; p0++) {
        // Visit each orbit once: from its first member in memory order, marking
        // the others as seen.
        const flat = p0 + s1 * p1 + s2 * p2;
        if (seen !== null && seen[flat >>> 5] & (1 << (flat & 31))) continue;
        // Pool the orbit's measured members. Each distinct member appears |stabilizer| times.
        let stab = 0, total = 0, valid = 0, errors = 0;
        for (let g = 0, o = 0; g < Ng; g++, o += 9) {
          const q0 = M[o] * p0 + M[o + 1] * p1 + M[o + 2] * p2 + T[3 * g];
          const q1 = M[o + 3] * p0 + M[o + 4] * p1 + M[o + 5] * p2 + T[3 * g + 1];
          const q2 = M[o + 6] * p0 + M[o + 7] * p1 + M[o + 8] * p2 + T[3 * g + 2];
          img[3 * g] = q0; img[3 * g + 1] = q1; img[3 * g + 2] = q2;
          if (q0 === p0 && q1 === p1 && q2 === p2) stab++;
          if (q0 < 0 || q1 < 0 || q2 < 0 || q0 >= n0 || q1 >= n1 || q2 >= n2) continue;
          const f = q0 + s1 * q1 + s2 * q2, v = volume[f];
          if (seen !== null) seen[f >>> 5] |= 1 << (f & 31);
          if (v === v && (mask === null || !mask[f])) {
            total += v;
            valid++;
            if (variance !== null) errors += variance[f];
          }
        }
        if (!valid) continue;
        orbits++;
        const mean = total / valid;
        // Variance of the mean over the distinct measured members.
        const error = variance !== null ? errors * stab / (valid * valid) : 0;
        // Sub-cells hit by the orbit's distinct members, per shell. Neighbouring
        // sub-cells mostly share a shell: |Q|² is first compared with its bounds.
        let m = 0, k = -1, lo = Infinity, hi = -Infinity;
        for (let g = 0; g < Ng; g++) {
          const x = c0[0] + img[3 * g] * w[0], y = c0[1] + img[3 * g + 1] * w[1], z = c0[2] + img[3 * g + 2] * w[2];
          // |Q|² at x + δ: xᵀGx + 2 δᵀGx + δᵀGδ.
          const gx = G[0] * x + G[1] * y + G[2] * z, gy = G[3] * x + G[4] * y + G[5] * z, gz = G[6] * x + G[7] * y + G[8] * z;
          const xx = x * gx + y * gy + z * gz;
          for (let a = 0, j = 0; a < split; a++) {
            const ta = xx + ex[a] * gx;
            for (let c1 = 0; c1 < split; c1++) {
              const tb = ta + ey[c1] * gy;
              for (let c2 = 0; c2 < split; c2++, j++) {
                const q2 = tb + ez[c2] * gz + dd[j];
                if (!(q2 >= lo && q2 < hi)) {
                  // sqrt of a rounding-negative |Q|² near 0 is NaN, and NaN | 0 is shell 0.
                  const b = (Math.sqrt(q2) * inv) | 0;
                  if (b >= count) { lo = Infinity; hi = -Infinity; continue; }
                  lo = b === 0 ? -Infinity : bound[b];
                  hi = bound[b + 1];
                  k = 0;
                  while (k < m && shells[k] !== b) k++;
                  if (k === m) { shells[m] = b; hits[m] = 0; m++; }
                }
                hits[k]++;
              }
            }
          }
        }
        // Every distinct member, split into S sub-cells, has weight 1 in all.
        const part = 1 / (stab * S);
        for (let k = 0; k < m; k++) {
          const b = shells[k], c = hits[k] * part;
          weight[b] += c;
          sum[b] += c * mean;
          spread[b] += c * c * error;
        }
      }
    }
    onProgress?.((p2 + 1) / n2);
  }

  const intensity = new Float64Array(count), sigma = new Float64Array(count), coverage = new Float64Array(count);
  for (let b = 0; b < count; b++) {
    const shell = (4 * Math.PI / 3) * ((b + 1) ** 3 - b ** 3) * dq ** 3;
    intensity[b] = weight[b] > 0 ? sum[b] / weight[b] : NaN;
    sigma[b] = weight[b] > 0 && variance !== null ? Math.sqrt(spread[b]) / weight[b] : NaN;
    coverage[b] = weight[b] * plan.voxel / shell;
  }
  const edges = Float64Array.from({ length: count + 1 }, (_, b) => b * dq);
  return { edges, intensity, sigma, voxels: weight, coverage, orbits };
}

function det3(m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
}
