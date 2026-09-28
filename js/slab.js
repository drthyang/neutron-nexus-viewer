// Slab averaging with symmetry pooling on an in-memory volume.
//
// Generalizes average_slab() from the original Python viewer: every output
// pixel is the equal-weight mean of the unique finite, unmasked source voxels
// in the union of the symmetry orbits of the slab voxels under that pixel.
// With the group {1, -1} this is exactly the original inversion averaging.
// Masked and non-finite voxels are stored as NaN, so "valid" is `v === v`.
//
// Conventions: `shape` is the storage (C-order) shape and display dimension d
// lives on storage axis 2 - d (Mantid writes signal with axes "D2:D1:D0").
// Operations are integer affine maps on display-dimension bin indices,
// i' = M i + t, as produced by indexMaps() in symmetry.js.

export const IDENTITY_MAP = { M: Int32Array.of(1, 0, 0, 0, 1, 0, 0, 0, 1), t: Int32Array.of(0, 0, 0) };

/** Bin indices whose centers lie within thickness/2 of center. */
export function selectBins(edges, center, thickness) {
  const ids = [];
  for (let i = 0; i + 1 < edges.length; i++) {
    if (Math.abs((edges[i] + edges[i + 1]) / 2 - center) < thickness / 2 + 1e-5) ids.push(i);
  }
  return ids;
}

/**
 * Average the slab `ids` along display dimension `fixed`, pooling each slab
 * voxel's symmetry orbit. Returns row-major values (rows follow the larger
 * remaining dimension), per-pixel unique voxel counts, and the filled fraction.
 *
 * Orbits are either identical or disjoint, so the pooled set is the disjoint
 * union of the distinct orbits met in the slab column; an orbit is identified
 * by its smallest in-range flat index. Within an orbit, every distinct image
 * appears |stabilizer| times among the group images, so sums over all images
 * are divided by the stabilizer order.
 */
export function averageSlab(volume, shape, fixed, ids, maps = [IDENTITY_MAP]) {
  const n = [shape[2], shape[1], shape[0]], strides = [1, n[0], n[0] * n[1]];
  const [x, y] = [0, 1, 2].filter((d) => d !== fixed);
  const R = n[y], C = n[x], S = ids.length;
  const values = new Float32Array(R * C), counts = new Uint32Array(R * C);
  // Per-column accumulators for the current row (and slab plane).
  const total = new Float64Array(C), count = new Float64Array(C), seen = new Int32Array(C);
  const sum = new Float64Array(C), valid = new Int32Array(C), stab = new Int32Array(C), rep = new Int32Array(C);
  const reps = new Int32Array(C * S);
  const ex = [0, 0, 0];
  ex[x] = 1;
  let filled = 0;

  for (let row = 0; row < R; row++) {
    total.fill(0); count.fill(0); seen.fill(0);
    for (let s = 0; s < S; s++) {
      sum.fill(0); valid.fill(0); stab.fill(0); rep.fill(0x7fffffff);
      const p = [0, 0, 0];
      p[y] = row;
      p[fixed] = ids[s];
      // Along a row only the x index changes, so each image q = M p + t moves
      // by column x of M; stream through the row one operation at a time.
      for (const { M, t } of maps) {
        const q = [0, 1, 2].map((i) => M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + t[i]);
        const d = [M[x], M[3 + x], M[6 + x]];
        const [fixLo, fixHi] = solveZero(q.map((qi, i) => qi - p[i]), d.map((di, i) => di - ex[i]), C);
        for (let col = fixLo; col < fixHi; col++) stab[col]++;
        const [lo, hi] = inRange(q, d, n, C);
        const df = d[0] + d[1] * strides[1] + d[2] * strides[2];
        let f = q[0] + q[1] * strides[1] + q[2] * strides[2] + lo * df;
        for (let col = lo; col < hi; col++, f += df) {
          if (f < rep[col]) rep[col] = f;
          const v = volume[f];
          if (v === v) { sum[col] += v; valid[col]++; }
        }
      }
      for (let col = 0; col < C; col++) {
        const r = rep[col], base = col * S;
        let known = false;
        for (let k = 0; k < seen[col]; k++) if (reps[base + k] === r) { known = true; break; }
        if (known) continue;
        reps[base + seen[col]++] = r;
        total[col] += sum[col] / stab[col];
        count[col] += valid[col] / stab[col];
      }
    }
    for (let col = 0, out = row * C; col < C; col++, out++) {
      counts[out] = count[col];
      if (count[col]) { values[out] = total[col] / count[col]; filled++; } else values[out] = NaN;
    }
  }
  return { values, counts, rows: R, cols: C, coverage: filled / values.length };
}

/** Columns c in [0, C) where e + d*c is zero in every component, as [lo, hi). */
function solveZero(e, d, C) {
  let lo = 0, hi = C;
  for (let i = 0; i < 3; i++) {
    if (d[i] === 0) {
      if (e[i] !== 0) return [0, 0];
    } else {
      const c = -e[i] / d[i];
      if (!Number.isInteger(c)) return [0, 0];
      lo = Math.max(lo, c);
      hi = Math.min(hi, c + 1);
    }
  }
  return lo < hi ? [lo, hi] : [0, 0];
}

/** Columns c in [0, C) where 0 <= q + d*c < n in every component, as [lo, hi). */
function inRange(q, d, n, C) {
  let lo = 0, hi = C;
  for (let i = 0; i < 3; i++) {
    if (d[i] === 0) {
      if (q[i] < 0 || q[i] >= n[i]) return [0, 0];
    } else if (d[i] > 0) {
      lo = Math.max(lo, Math.ceil(-q[i] / d[i]));
      hi = Math.min(hi, Math.floor((n[i] - 1 - q[i]) / d[i]) + 1);
    } else {
      lo = Math.max(lo, Math.ceil((n[i] - 1 - q[i]) / d[i]));
      hi = Math.min(hi, Math.floor(-q[i] / d[i]) + 1);
    }
  }
  return lo < hi ? [lo, hi] : [0, 0];
}
