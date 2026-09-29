// Line cuts: the 1-D profile of the volume along a line, averaged over a rod.
//
// A cut runs from a to b in display coordinates (it is drawn in a slice, so it
// lies in that slice's plane). It is binned along its dominant axis `dom`, the
// display axis it changes most along, so its points are values of that
// coordinate: bin k is centred at a[dom] + k step. A voxel belongs to bin k when
// its centre lies within `radius` of the line, inside a rod around it, and
// projects onto the line within step/2 of the bin centre. Distances are measured
// in reciprocal space, x -> T x with T the display->Cartesian matrix (row-major,
// see cartesianBasis()), so the rod is round across the line whatever the
// lattice, and its width is the same in every direction across it.
//
// As in the slices, each bin is the equal-weight mean of the unique finite,
// unmasked voxels in the union of the symmetry orbits of the voxels in it. With
// per-voxel variances, σ is that of a mean of independent voxels, √(Σσ²)/N.
//
// Conventions as in slab.js: `shape` is the storage (C-order) shape, display
// dimension d lives on storage axis 2 - d, and symmetry maps are integer affine
// maps i' = M i + t on display-dimension bin indices (from indexMaps()).

import { IDENTITY_MAP } from './slab.js';

/** The axis a cut from a to b changes most along (the first of equals), for points of any dimension. */
export const cutAxis = (a, b) => a.reduce((best, _, d) => (Math.abs(b[d] - a[d]) > Math.abs(b[best] - a[best]) ? d : best), 0);

/**
 * Bin edges of a cut along its dominant axis: bins of width `step` centred on
 * a[dom], a[dom] + step, … up to the one nearest b[dom], so a cut between two
 * voxel centres samples each voxel column once. `a[dom]` must not exceed `b[dom]`.
 */
export function cutEdges(a, b, dom, step, maxBins = 20000) {
  const count = Math.round((b[dom] - a[dom]) / step) + 1;
  if (!(count <= maxBins)) throw new Error(`A cut of ${count} points is too fine: increase the step.`);
  return Float64Array.from({ length: count + 1 }, (_, k) => a[dom] + (k - 0.5) * step);
}

const apply = (T, x) => [0, 1, 2].map((i) => T[3 * i] * x[0] + T[3 * i + 1] * x[1] + T[3 * i + 2] * x[2]);
const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
const unit = (v) => { const size = Math.hypot(...v); return v.map((x) => x / size); };
/** The unit normal, in reciprocal space, of the planes where display coordinate `fixed` is constant (Ti = T⁻¹). */
const normalOf = (Ti, fixed) => unit([Ti[3 * fixed], Ti[3 * fixed + 1], Ti[3 * fixed + 2]]);

/** Inverse of a row-major 3×3 matrix. */
function inverse(T) {
  const [a, b, c, d, e, f, g, h, i] = T;
  const A = e * i - f * h, B = f * g - d * i, C = d * h - e * g, det = a * A + b * B + c * C;
  return [A, c * h - b * i, b * f - c * e, B, a * i - c * g, c * d - a * f, C, b * g - a * h, a * e - b * d].map((v) => v / det);
}

/**
 * The axis of a cut's rod: its start A and unit direction in reciprocal space,
 * the change of its axis coordinate per unit length (`rate`), and the display
 * point at a value s of that coordinate.
 */
function rodAxis({ a, b, dom, T }) {
  const A = apply(T, a), B = apply(T, b), D = B.map((v, i) => v - A[i]), length = Math.hypot(...D);
  if (!(length > 0)) throw new Error('The cut has no length.');
  const rate = (b[dom] - a[dom]) / length;
  if (!(rate > 0)) throw new Error('A cut must run toward larger values of its axis.');
  const at = (s) => a.map((ad, d) => ad + ((s - a[dom]) / (b[dom] - a[dom])) * (b[d] - ad));
  return { A, dir: D.map((v) => v / length), rate, at };
}

/**
 * Where a cut's rod meets the plane through its line normal to display axis
 * `fixed` (its slice): the corners of that band, along the line from its first
 * bin edge to its last and `radius` to either side, in display coordinates and
 * drawing order.
 */
export function cutBand(cut, fixed) {
  const { a, dom, edges, radius, T } = cut, { A, dir, rate } = rodAxis(cut), Ti = inverse(T);
  // Across the line within the plane: normal to the line and to the plane's normal.
  const across = unit(cross(normalOf(Ti, fixed), dir));
  const t0 = (edges[0] - a[dom]) / rate, t1 = (edges[edges.length - 1] - a[dom]) / rate;
  return [[t0, -radius], [t1, -radius], [t1, radius], [t0, radius]]
    .map(([t, side]) => apply(Ti, A.map((v, i) => v + t * dir[i] + side * across[i])));
}

/**
 * Where a cut's rod crosses the plane on which display coordinate `fixed` is
 * `value` (a slice the cut does not lie in): `at`, the point on its line, and
 * `ring`, the outline of the rod's section there, an ellipse `radius` across
 * the line and radius / cos θ along it (θ between the line and the plane's
 * normal), both in display coordinates. Null when the line runs along the
 * plane or meets it beyond the rod's ends.
 */
export function cutCrossing(cut, fixed, value, points = 48) {
  const { a, b, dom, edges, radius, T } = cut;
  if (!(Math.abs(b[fixed] - a[fixed]) > 1e-12)) return null;
  const f = (value - a[fixed]) / (b[fixed] - a[fixed]), at = a.map((ad, d) => ad + f * (b[d] - ad));
  if (!(at[dom] >= edges[0] && at[dom] <= edges[edges.length - 1])) return null;
  const { dir } = rodAxis(cut), Ti = inverse(T), normal = normalOf(Ti, fixed);
  const cos = dir[0] * normal[0] + dir[1] * normal[1] + dir[2] * normal[2];
  // In the plane: u along the line's shadow on it (any direction when the line crosses square on), v across.
  const shadow = dir.map((v, i) => v - cos * normal[i]);
  const u = Math.hypot(...shadow) > 1e-9 ? unit(shadow) : unit(cross(normal, Math.abs(normal[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
  const v = cross(normal, u), X = apply(T, at), long = radius / Math.abs(cos);
  const ring = Array.from({ length: points }, (_, k) => {
    const c = long * Math.cos((2 * Math.PI * k) / points), s = radius * Math.sin((2 * Math.PI * k) / points);
    return apply(Ti, X.map((x, i) => x + c * u[i] + s * v[i]));
  });
  return { at, ring };
}

/**
 * The line cut `cut` = { a, b, dom, edges, radius, T } of `volume`, whose axes
 * are `dims`. Returns per bin the intensity, σ (NaN without variances) and the
 * number of distinct voxels with data.
 */
export function lineCut(volume, shape, dims, cut, maps = [IDENTITY_MAP], mask = null, invert = false, variance = null) {
  const { a, dom, edges, radius, T } = cut, { A, dir, rate, at } = rodAxis(cut), Ti = inverse(T);
  const n = [shape[2], shape[1], shape[0]], strides = [1, n[0], n[0] * n[1]];
  const count = edges.length - 1, step = (edges[count] - edges[0]) / count;
  // The voxels to test: the box around the rod's axis, from its first bin edge
  // to its last, widened on each display axis by how far the radius reaches
  // along it (the radius times that coordinate's gradient), and one bin more.
  const [p0, p1] = [at(edges[0]), at(edges[count])];
  const axes = [0, 1, 2].map((d) => {
    const e = dims[d].edges, w = (e[e.length - 1] - e[0]) / (e.length - 1);
    const reach = radius * Math.hypot(Ti[3 * d], Ti[3 * d + 1], Ti[3 * d + 2]);
    return {
      e0: e[0], w,
      lo: Math.max(0, Math.floor((Math.min(p0[d], p1[d]) - reach - e[0]) / w) - 1),
      hi: Math.min(n[d] - 1, Math.floor((Math.max(p0[d], p1[d]) + reach - e[0]) / w) + 1),
    };
  });

  const Ng = maps.length, img = [0, 0, 0], p = [0, 0, 0], x = [0, 0, 0];
  const total = new Float64Array(count), valid = new Float64Array(count), spread = new Float64Array(count);
  const seen = Array.from({ length: count }, () => new Set());
  const reach2 = (radius * (1 + 1e-9)) ** 2 + 1e-18;
  for (p[2] = axes[2].lo; p[2] <= axes[2].hi; p[2]++) {
    for (p[1] = axes[1].lo; p[1] <= axes[1].hi; p[1]++) {
      for (p[0] = axes[0].lo; p[0] <= axes[0].hi; p[0]++) {
        for (let d = 0; d < 3; d++) x[d] = axes[d].e0 + (p[d] + 0.5) * axes[d].w;
        const X = apply(T, x), r = X.map((v, i) => v - A[i]);
        const t = r[0] * dir[0] + r[1] * dir[1] + r[2] * dir[2];
        if (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] - t * t > reach2) continue;
        const k = Math.floor((a[dom] + t * rate - edges[0]) / step + 1e-9);
        if (k < 0 || k >= count) continue;
        // The voxel's orbit: its smallest in-range flat index names it, and every
        // distinct member appears |stabilizer| times among the group images.
        let rep = Infinity, stab = 0, sum = 0, n0 = 0, var0 = 0;
        for (const { M, t: shift } of maps) {
          for (let i = 0; i < 3; i++) img[i] = M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + shift[i];
          if (img[0] === p[0] && img[1] === p[1] && img[2] === p[2]) stab++;
          if (img[0] < 0 || img[1] < 0 || img[2] < 0 || img[0] >= n[0] || img[1] >= n[1] || img[2] >= n[2]) continue;
          const f = img[0] + img[1] * strides[1] + img[2] * strides[2];
          if (f < rep) rep = f;
          const value = volume[f];
          if (value === value && (mask === null || (mask[f] !== 0) === invert)) {
            sum += value;
            n0++;
            if (variance) var0 += variance[f];
          }
        }
        if (Ng > 1) {
          if (seen[k].has(rep)) continue;
          seen[k].add(rep);
        }
        total[k] += sum / stab;
        valid[k] += n0 / stab;
        spread[k] += var0 / stab;
      }
    }
  }
  const intensity = new Float64Array(count), sigma = new Float64Array(count), voxels = new Float64Array(count);
  for (let k = 0; k < count; k++) {
    voxels[k] = valid[k];
    intensity[k] = valid[k] ? total[k] / valid[k] : NaN;
    sigma[k] = valid[k] && variance ? Math.sqrt(spread[k]) / valid[k] : NaN;
  }
  return { edges: Float64Array.from(edges), intensity, sigma, voxels };
}
