// User masks on the unsymmetrized volume, applied before symmetry averaging.
//
// Reduced volumes often carry spurious high values where detector edges (and
// their weak normalization) land in reciprocal space: along the boundaries of
// the measured region. Two complementary masks target them:
//
//   EDGE     measured voxels within `radius` voxels (box distance) of an
//            unmeasured voxel, i.e. an erosion of the coverage;
//   OUTLIER  voxels that exceed the median of their symmetry-equivalent voxels
//            by more than k standard deviations: the robust spread 1.4826 * MAD,
//            but at least the orbit's counting error, so that counting noise in
//            low-count data (where MAD is often 0) is not taken for outliers.
//
// The mask is a Uint8Array in storage order with these bits; volume values are
// never modified, so masks can be changed or cleared at any time.

export const EDGE = 1;
export const OUTLIER = 2;

/** Flag measured voxels within `radius` voxels of an unmeasured (NaN) voxel. */
export function edgeMask(volume, shape, radius, mask = new Uint8Array(volume.length)) {
  if (radius < 1) return mask;
  // Dilate the unmeasured set with a (2r+1)^3 box, one axis at a time.
  const near = new Uint8Array(volume.length);
  for (let i = 0; i < volume.length; i++) near[i] = volume[i] !== volume[i] ? 1 : 0;
  const strides = [shape[1] * shape[2], shape[2], 1];
  for (let axis = 0; axis < 3; axis++) {
    const L = shape[axis], s = strides[axis], line = new Uint8Array(L);
    const others = [0, 1, 2].filter((a) => a !== axis);
    for (let a = 0; a < shape[others[0]]; a++) {
      for (let b = 0; b < shape[others[1]]; b++) {
        const base = a * strides[others[0]] + b * strides[others[1]];
        for (let i = 0; i < L; i++) line[i] = near[base + i * s];
        let count = 0;
        for (let i = 0; i < Math.min(radius, L); i++) count += line[i];
        for (let i = 0; i < L; i++) {
          if (i + radius < L) count += line[i + radius];
          if (i - radius - 1 >= 0) count -= line[i - radius - 1];
          near[base + i * s] = count > 0 ? 1 : 0;
        }
      }
    }
  }
  for (let i = 0; i < volume.length; i++) if (near[i] && volume[i] === volume[i]) mask[i] |= EDGE;
  return mask;
}

/**
 * Flag voxels more than k σ above the median of their orbit's valid voxels
 * (those not NaN and not already masked). Orbits with fewer than `minMembers`
 * valid voxels are left alone. `maps` are fine-grid index maps, identity first.
 *
 * σ is the larger of the robust spread 1.4826·MAD and the orbit's typical
 * counting error √median(σᵥ²), with σᵥ² from `variance` (per voxel, storage
 * order, NaN unknown) or, without it, the counts themselves when every member
 * is a non-negative integer (Poisson). When that is 0 too (more than half the
 * members equal, as zero-count voxels are), a voxel is tested against its own
 * σᵥ instead, so one bright voxel among empty ones is still caught; without
 * any error estimate such an orbit is skipped. Without the floor, σ = 0 flagged
 * every voxel above the median: a quarter of a pure Poisson volume at 0.3 counts.
 */
export function outlierMask(volume, shape, maps, k, mask = new Uint8Array(volume.length), { minMembers = 3, variance = null, onProgress = () => {} } = {}) {
  const n0 = shape[2], n1 = shape[1], n2 = shape[0], s1 = n0, s2 = n0 * n1;
  const G = maps.length, flats = new Int32Array(G), values = new Float64Array(G), errors = new Float64Array(G), scratch = new Float64Array(G);
  const usable = (f) => volume[f] === volume[f] && !mask[f];
  const flagged = [];
  for (let i2 = 0, p = 0; i2 < n2; i2++) {
    for (let i1 = 0; i1 < n1; i1++) {
      for (let i0 = 0; i0 < n0; i0++, p++) {
        if (!usable(p)) continue;
        // Visit each orbit once, from its smallest usable flat index.
        let m = 0, rep = true;
        for (let g = 0; g < G && rep; g++) {
          const { M, t } = maps[g];
          const q0 = M[0] * i0 + M[1] * i1 + M[2] * i2 + t[0];
          const q1 = M[3] * i0 + M[4] * i1 + M[5] * i2 + t[1];
          const q2 = M[6] * i0 + M[7] * i1 + M[8] * i2 + t[2];
          if (q0 < 0 || q1 < 0 || q2 < 0 || q0 >= n0 || q1 >= n1 || q2 >= n2) continue;
          const f = q0 + q1 * s1 + q2 * s2;
          if (!usable(f)) continue;
          if (f < p) rep = false;
          else if (!flats.subarray(0, m).includes(f)) flats[m++] = f;
        }
        if (!rep || m < minMembers) continue;
        let counts = !variance;
        for (let j = 0; j < m; j++) {
          const v = volume[flats[j]];
          scratch[j] = values[j] = v;
          if (counts && !(v >= 0 && Number.isInteger(v))) counts = false;
        }
        const median = medianOf(scratch.subarray(0, m));
        for (let j = 0; j < m; j++) scratch[j] = Math.abs(values[j] - median);
        let sigma = 1.4826 * medianOf(scratch.subarray(0, m));
        // Counting errors: σᵥ² from the file, or the count itself; NaN when unknown.
        let known = 0;
        for (let j = 0; j < m; j++) {
          const v2 = variance ? variance[flats[j]] : counts ? values[j] : NaN;
          errors[j] = v2 >= 0 ? Math.sqrt(v2) : NaN;
          if (v2 >= 0) scratch[known++] = v2;
        }
        if (known) sigma = Math.max(sigma, Math.sqrt(medianOf(scratch.subarray(0, known))));
        for (let j = 0; j < m; j++) {
          const s = sigma > 0 ? sigma : errors[j];
          if (s > 0 && values[j] - median > k * s) flagged.push(flats[j]);
        }
      }
    }
    if (i2 % 16 === 0) onProgress((i2 + 1) / n2);
  }
  for (const f of flagged) mask[f] |= OUTLIER;
  onProgress(1);
  return mask;
}

/** Median of a typed array; sorts it in place. */
function medianOf(a) {
  a.sort();
  const n = a.length;
  return n % 2 ? a[(n - 1) / 2] : (a[n / 2 - 1] + a[n / 2]) / 2;
}

/** Counts of measured voxels per mask bit. */
export function maskStats(volume, mask) {
  let measured = 0, edge = 0, outlier = 0;
  for (let i = 0; i < volume.length; i++) {
    if (volume[i] !== volume[i]) continue;
    measured++;
    if (mask[i] & EDGE) edge++;
    else if (mask[i] & OUTLIER) outlier++;
  }
  return { measured, edge, outlier };
}
