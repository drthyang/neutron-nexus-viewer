// Line cuts: the 1-D profile of the volume along a line in a slice plane.
//
// A cut runs from a to b in the plane of the slice through display dims x and y,
// within the slab of that slice (the bins of dim `fixed` whose centres lie within
// thickness/2 of `center`). It is binned along its dominant axis `dom` (0: the
// plane's x axis, 1: its y axis), the one it changes more along, so its points
// are values of that coordinate: bin k is centred at a[dom] + k step. A voxel
// belongs to the bin its centre projects into, if the centre lies within `half`
// of the line, measured across it in the plane's Cartesian metric
// (planeGeometry(): lengths lx and ly per unit coordinate, cos of their angle).
//
// As in the slices, each bin is the equal-weight mean of the unique finite,
// unmasked voxels in the union of the symmetry orbits of the voxels in it. With
// per-voxel variances, σ is that of a mean of independent voxels, √(Σσ²)/N.
//
// Conventions as in slab.js: `shape` is the storage (C-order) shape, display
// dimension d lives on storage axis 2 - d, and symmetry maps are integer affine
// maps i' = M i + t on display-dimension bin indices (from indexMaps()).

import { IDENTITY_MAP, selectBins } from './slab.js';

/** The axis a cut from a to b (in-plane coordinates) is binned along: 0 (x) or 1 (y). */
export const cutAxis = (a, b) => (Math.abs(b[1] - a[1]) > Math.abs(b[0] - a[0]) ? 1 : 0);

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

/**
 * The plane's Cartesian metric: `cart` maps plane coordinates (u, v) to lengths,
 * `plane` maps back. Lengths per unit coordinate lx and ly, cos of their angle.
 */
function metric({ lx, ly, cos }) {
  const sin = Math.sqrt(Math.max(0, 1 - cos * cos));
  return {
    cart: (u, v) => [lx * u + ly * cos * v, ly * sin * v],
    plane: (X, Y) => { const v = Y / (ly * sin); return [(X - ly * cos * v) / lx, v]; },
  };
}

/**
 * The corners of a cut's band in plane coordinates: along the line from its
 * first bin edge to its last, `half` to either side, in drawing order.
 */
export function cutBand({ a, b, dom, edges, half, geometry }) {
  const { cart, plane } = metric(geometry);
  const [ax, ay] = cart(a[0], a[1]), [bx, by] = cart(b[0], b[1]);
  const length = Math.hypot(bx - ax, by - ay), tx = (bx - ax) / length, ty = (by - ay) / length;
  const rate = (b[dom] - a[dom]) / length, t0 = (edges[0] - a[dom]) / rate, t1 = (edges[edges.length - 1] - a[dom]) / rate;
  return [[t0, -half], [t1, -half], [t1, half], [t0, half]].map(([t, side]) => plane(ax + t * tx - side * ty, ay + t * ty + side * tx));
}

/**
 * The line cut `cut` = { fixed, x, y, center, thickness, a, b, dom, edges, half,
 * geometry: { lx, ly, cos } } of `volume`, whose axes are `dims`. Returns per
 * bin the intensity, σ (NaN without variances) and the number of distinct voxels
 * with data, and the number of slab bins used.
 */
export function lineCut(volume, shape, dims, cut, maps = [IDENTITY_MAP], mask = null, invert = false, variance = null) {
  const { fixed, x, y, center, thickness, a, b, dom, edges, half } = cut;
  const ids = selectBins(dims[fixed].edges, center, thickness);
  if (!ids.length) throw new Error('No bins selected: increase the thickness or move the center.');
  const n = [shape[2], shape[1], shape[0]], strides = [1, n[0], n[0] * n[1]];
  const { cart } = metric(cut.geometry);
  const [ax, ay] = cart(a[0], a[1]), [bx, by] = cart(b[0], b[1]);
  const length = Math.hypot(bx - ax, by - ay);
  if (!(length > 0)) throw new Error('The cut has no length.');
  const tx = (bx - ax) / length, ty = (by - ay) / length;
  const count = edges.length - 1, step = (edges[count] - edges[0]) / count, rate = (b[dom] - a[dom]) / length;
  if (!(rate > 0)) throw new Error('A cut must run toward larger values of its axis.');

  // Bin centres along each in-plane axis (uniform bins), and the range of each
  // that the band covers: its corners, one bin wider.
  const axis = (d) => {
    const e = dims[d].edges, w = (e[e.length - 1] - e[0]) / (e.length - 1);
    return { e0: e[0], w, count: e.length - 1 };
  };
  const U = axis(x), V = axis(y), corners = cutBand(cut);
  const span = ({ e0, w, count: m }, values) => [
    Math.max(0, Math.floor((Math.min(...values) - e0) / w) - 1),
    Math.min(m - 1, Math.floor((Math.max(...values) - e0) / w) + 1),
  ];
  const [i0, i1] = span(U, corners.map((c) => c[0])), [j0, j1] = span(V, corners.map((c) => c[1]));

  const Ng = maps.length, img = [0, 0, 0], p = [0, 0, 0];
  const total = new Float64Array(count), valid = new Float64Array(count), spread = new Float64Array(count);
  const seen = Array.from({ length: count }, () => new Set());
  const tolerance = 1e-9 * Math.max(1, half);
  for (const k0 of ids) {
    p[fixed] = k0;
    for (let j = j0; j <= j1; j++) {
      p[y] = j;
      const v = V.e0 + (j + 0.5) * V.w;
      for (let i = i0; i <= i1; i++) {
        const u = U.e0 + (i + 0.5) * U.w;
        const [X, Y] = cart(u, v), rx = X - ax, ry = Y - ay;
        if (Math.abs(ry * tx - rx * ty) > half + tolerance) continue;
        const k = Math.floor((a[dom] + (rx * tx + ry * ty) * rate - edges[0]) / step + 1e-9);
        if (k < 0 || k >= count) continue;
        p[x] = i;
        // The voxel's orbit: its smallest in-range flat index names it, and every
        // distinct member appears |stabilizer| times among the group images.
        let rep = Infinity, stab = 0, sum = 0, n0 = 0, var0 = 0;
        for (const { M, t } of maps) {
          for (let r = 0; r < 3; r++) img[r] = M[3 * r] * p[0] + M[3 * r + 1] * p[1] + M[3 * r + 2] * p[2] + t[r];
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
  return { edges: Float64Array.from(edges), intensity, sigma, voxels, bins: ids.length, slab: [dims[fixed].edges[ids[0]], dims[fixed].edges[ids.at(-1) + 1]] };
}
