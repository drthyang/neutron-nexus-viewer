// Reduction of Rigaku Oxford Diffraction frames to an HKL volume, in three passes
// over the frames (each frame is read and decoded once per pass, so memory stays at
// a few detector-sized arrays plus the output grid):
//   1. per-pixel sums -> static detector mask (border, chip boundaries, beamstop
//      umbra and penumbra, dead pixels) and a per-run background
//   2. 3-D peak search (strong pixels joined across neighbouring pixels and frames)
//   3. after refining the geometry against the peaks (rigaku-geometry.js), every
//      unmasked pixel of every frame is binned into the grid
// With a WorkerPool (rigaku-pool.js) the per-frame work of each pass runs in parallel
// workers (rigaku-map-worker.js) and is merged so that the result is bit-identical to the
// serial path: integer sums are exact in any order, the peak search runs one whole run per
// worker, and the gridding contributions are added in frame and pixel order by one thread.
// The method, its parameters and its validation against the reference Python
// implementation are described in docs/METHOD.md (Rigaku reduction).

import {
  CELL_SYSTEMS, DEG, M_CRYSALIS, bMatrix, cellFromUB, crystalSystem, det3, gonio, headerGeometry, inv3, levenbergMarquardt,
  mul, mulv, pixelToLab, polarU, predict, prepare, samplePosition, smallRot, transformCell, transpose,
} from './rigaku-geometry.js';

// ---- Small image tools ----------------------------------------------------------

/** k-th smallest of arr[0..n) (quickselect, in place). */
function select(arr, n, k) {
  let lo = 0, hi = n - 1;
  while (hi > lo) {
    const pivot = arr[(lo + hi) >> 1];
    let i = lo, j = hi;
    while (i <= j) {
      while (arr[i] < pivot) i++;
      while (arr[j] > pivot) j--;
      if (i <= j) { const t = arr[i]; arr[i] = arr[j]; arr[j] = t; i++; j--; }
    }
    if (k <= j) hi = j; else if (k >= i) lo = i; else return arr[k];
  }
  return arr[k];
}

/**
 * Smooth quantile filter: the q-quantile of the finite values in a `size` x `size`
 * window, evaluated every `stride` pixels and interpolated bilinearly. Values where
 * `skip` is set are ignored.
 */
export function quantileFilter(img, nx, ny, size, q, skip = null, stride = 4) {
  const gx = Math.ceil(nx / stride) + 1, gy = Math.ceil(ny / stride) + 1, h = size >> 1;
  const coarse = new Float64Array(gx * gy), buf = new Float64Array(size * size);
  for (let b = 0; b < gy; b++) {
    const yc = Math.min(ny - 1, b * stride);
    for (let a = 0; a < gx; a++) {
      const xc = Math.min(nx - 1, a * stride);
      let n = 0;
      for (let y = Math.max(0, yc - h); y <= Math.min(ny - 1, yc + h); y++) {
        for (let x = Math.max(0, xc - h); x <= Math.min(nx - 1, xc + h); x++) {
          const k = y * nx + x;
          if (skip && skip[k]) continue;
          const v = img[k];
          if (v === v) buf[n++] = v;
        }
      }
      coarse[b * gx + a] = n ? select(buf, n, Math.min(n - 1, Math.floor(q * (n - 1) + 0.5))) : NaN;
    }
  }
  const out = new Float64Array(nx * ny);
  for (let y = 0; y < ny; y++) {
    const fy = Math.min(y / stride, gy - 1), b0 = Math.min(Math.floor(fy), gy - 2), ty = fy - b0;
    for (let x = 0; x < nx; x++) {
      const fx = Math.min(x / stride, gx - 1), a0 = Math.min(Math.floor(fx), gx - 2), tx = fx - a0;
      const v00 = coarse[b0 * gx + a0], v01 = coarse[b0 * gx + a0 + 1], v10 = coarse[(b0 + 1) * gx + a0], v11 = coarse[(b0 + 1) * gx + a0 + 1];
      out[y * nx + x] = (1 - ty) * ((1 - tx) * v00 + tx * v01) + ty * ((1 - tx) * v10 + tx * v11);
    }
  }
  return out;
}

/** Connected components (4-connectivity) of a boolean image; returns [labels, sizes]. */
export function label2D(on, nx, ny) {
  const lab = new Int32Array(nx * ny).fill(-1), sizes = [];
  const stack = new Int32Array(nx * ny);
  for (let s = 0; s < nx * ny; s++) {
    if (!on[s] || lab[s] >= 0) continue;
    const id = sizes.length;
    let top = 0, n = 0;
    stack[top++] = s;
    lab[s] = id;
    while (top) {
      const k = stack[--top];
      n++;
      const x = k % nx, y = (k - x) / nx;
      for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const u = x + dx, v = y + dy;
        if (u < 0 || v < 0 || u >= nx || v >= ny) continue;
        const q = v * nx + u;
        if (on[q] && lab[q] < 0) { lab[q] = id; stack[top++] = q; }
      }
    }
    sizes.push(n);
  }
  return [lab, sizes];
}

/** Binary dilation with the 4-neighbour cross, `iterations` times. */
export function dilate(on, nx, ny, iterations = 1) {
  let cur = Uint8Array.from(on);
  for (let it = 0; it < iterations; it++) {
    const next = Uint8Array.from(cur);
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const k = y * nx + x;
        if (cur[k]) continue;
        if ((x > 0 && cur[k - 1]) || (x < nx - 1 && cur[k + 1]) || (y > 0 && cur[k - nx]) || (y < ny - 1 && cur[k + nx])) next[k] = 1;
      }
    }
    cur = next;
  }
  return cur;
}

/** Euclidean distance (pixels) to the nearest set pixel (Felzenszwalb-Huttenlocher). */
export function distanceTransform(on, nx, ny) {
  const INF = 1e20, f = new Float64Array(Math.max(nx, ny)), d = new Float64Array(Math.max(nx, ny));
  const v = new Int32Array(Math.max(nx, ny)), z = new Float64Array(Math.max(nx, ny) + 1);
  const g = new Float64Array(nx * ny);
  const pass = (n) => {
    let k = 0;
    v[0] = 0; z[0] = -INF; z[1] = INF;
    for (let q = 1; q < n; q++) {
      let s;
      while (true) {
        s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        if (s <= z[k]) k--; else break;
      }
      k++; v[k] = q; z[k] = s; z[k + 1] = INF;
    }
    k = 0;
    for (let q = 0; q < n; q++) {
      while (z[k + 1] < q) k++;
      d[q] = (q - v[k]) ** 2 + f[v[k]];
    }
  };
  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) f[y] = on[y * nx + x] ? 0 : INF;
    pass(ny);
    for (let y = 0; y < ny; y++) g[y * nx + x] = d[y];
  }
  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) f[x] = g[y * nx + x];
    pass(nx);
    for (let x = 0; x < nx; x++) g[y * nx + x] = Math.sqrt(d[x]);
  }
  return g;
}

function median(arr) {
  const a = Float64Array.from(arr.filter((x) => x === x)).sort();
  return a.length ? (a.length % 2 ? a[(a.length - 1) / 2] : 0.5 * (a[a.length / 2 - 1] + a[a.length / 2])) : NaN;
}

// ---- Detector mask ----------------------------------------------------------------

/**
 * Static detector mask (1 = excluded) from the per-pixel sum over all frames and the
 * number of frames in which each pixel counted:
 *  - the outer 1-pixel border;
 *  - chip-boundary triplets: rows/columns whose median deviates by more than 8 % from
 *    a smooth trend (re-binned large boundary pixels), plus any `boundaries` given;
 *  - the beamstop umbra (sum < 20 % of the local median level, components of >= 50 px,
 *    dilated by 3) and penumbra (connected pixels below 75 % of the local 90th
 *    percentile, at most 26 px from the umbra, dilated by 2);
 *  - pixels that never counted.
 */
export function detectorMask(sum, counted, nx, ny, { boundaries = null } = {}) {
  const mask = new Uint8Array(nx * ny), info = {};
  for (let x = 0; x < nx; x++) mask[x] = mask[(ny - 1) * nx + x] = 1;
  for (let y = 0; y < ny; y++) mask[y * nx] = mask[y * nx + nx - 1] = 1;
  const level = quantileFilter(sum, nx, ny, 41, 0.5);
  const low = Uint8Array.from(sum, (v, k) => (v < 0.2 * level[k] ? 1 : 0));
  const [lab, sizes] = label2D(low, nx, ny);
  const umbraCore = Uint8Array.from(lab, (l) => (l >= 0 && sizes[l] >= 50 ? 1 : 0));
  let umbra = dilate(umbraCore, nx, ny, 3);
  info.umbraPixels = umbra.reduce((s, v) => s + v, 0);
  const ref = quantileFilter(sum, nx, ny, 61, 0.9, umbra);
  const pen = Uint8Array.from(sum, (v, k) => (v < 0.75 * ref[k] || umbra[k] ? 1 : 0));
  const [lab2] = label2D(pen, nx, ny);
  const keep = new Set();
  for (let k = 0; k < nx * ny; k++) if (umbra[k] && lab2[k] >= 0) keep.add(lab2[k]);
  const dist = distanceTransform(umbra, nx, ny);
  const pen2 = Uint8Array.from(lab2, (l, k) => (l >= 0 && keep.has(l) && dist[k] <= 26 ? 1 : 0));
  const shadow = dilate(pen2, nx, ny, 2);
  info.beamstopPixels = shadow.reduce((s, v) => s + v, 0);
  // chip boundaries from column/row medians outside the shadow
  const prof = (along) => {
    const n = along === 'col' ? nx : ny, m = along === 'col' ? ny : nx;
    const out = new Float64Array(n), tmp = [];
    for (let a = 0; a < n; a++) {
      tmp.length = 0;
      for (let b = 1; b < m - 1; b++) {
        const k = along === 'col' ? b * nx + a : a * nx + b;
        if (!shadow[k]) tmp.push(sum[k]);
      }
      out[a] = median(tmp);
    }
    return out;
  };
  const lines = (p) => {
    const n = p.length, bad = [];
    for (let a = 1; a < n - 1; a++) {
      const w = [];
      for (let b = Math.max(1, a - 7); b <= Math.min(n - 2, a + 7); b++) w.push(p[b]);
      const t = median(w);
      if (Math.abs(p[a] / t - 1) > 0.08) bad.push(a);
    }
    const out = new Set();
    let group = [];
    for (const i of [...bad, null]) {
      if (group.length && (i === null || i - group.at(-1) > 2)) {
        const c = Math.round(group.reduce((s, v) => s + v, 0) / group.length);
        out.add(c - 1); out.add(c); out.add(c + 1);
        group = [];
      }
      if (i !== null) group.push(i);
    }
    return [...out].filter((i) => i >= 0 && i < n).sort((a, b) => a - b);
  };
  let cols = lines(prof('col')), rows = lines(prof('row'));
  info.boundaryColumnsDetected = cols.slice();
  info.boundaryRowsDetected = rows.slice();
  if (boundaries) {
    cols = [...new Set([...cols, ...boundaries.columns])].sort((a, b) => a - b);
    rows = [...new Set([...rows, ...boundaries.rows])].sort((a, b) => a - b);
  }
  for (const c of cols) for (let y = 0; y < ny; y++) mask[y * nx + c] = 1;
  for (const r of rows) for (let x = 0; x < nx; x++) mask[r * nx + x] = 1;
  let dead = 0;
  for (let k = 0; k < nx * ny; k++) {
    if (shadow[k]) mask[k] = 1;
    if (!counted[k]) { mask[k] = 1; dead++; }
  }
  info.boundaryColumns = cols;
  info.boundaryRows = rows;
  info.deadPixels = dead;
  info.maskedPixels = mask.reduce((s, v) => s + v, 0);
  info.maskedFraction = info.maskedPixels / mask.length;
  return { mask, info };
}

// ---- Peak search ------------------------------------------------------------------

/**
 * Streaming 3-D peak search over the frames of one run. Strong pixels: counts >= 8 and
 * >= bg + 6 sqrt(bg + 1) (bg = the pixel's mean over the run), outside the mask; they are
 * joined across the 4 in-frame neighbours and the same pixel in the previous frame.
 * Components with >= 4 voxels and >= 150 net counts are kept; centroids are net-count
 * weighted, the frame coordinate at frame midpoints.
 */
export class PeakSearch {
  constructor(nx, ny, mask, bg) {
    Object.assign(this, { nx, ny, mask, bg });
    this.near = dilate(mask, nx, ny, 1);
    this.prev = new Int32Array(nx * ny).fill(-1);
    this.cur = new Int32Array(nx * ny).fill(-1);
    this.parent = [];
    this.st = []; // [w, wx, wy, wz, n, max, zmin, zmax, touches]
    this.z = 0;
  }

  find(a) {
    const p = this.parent;
    while (p[a] !== a) { p[a] = p[p[a]]; a = p[a]; }
    return a;
  }

  union(a, b) {
    a = this.find(a); b = this.find(b);
    if (a === b) return a;
    if (b < a) [a, b] = [b, a];
    this.parent[b] = a;
    const A = this.st[a], B = this.st[b];
    for (let i = 0; i < 5; i++) A[i] += B[i];
    A[5] = Math.max(A[5], B[5]); A[6] = Math.min(A[6], B[6]); A[7] = Math.max(A[7], B[7]); A[8] = A[8] || B[8];
    return a;
  }

  addFrame(img) {
    const { nx, ny, mask, bg, near, prev, cur } = this, z = this.z;
    cur.fill(-1);
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const k = y * nx + x;
        if (mask[k]) continue;
        const c = img[k], b = bg[k];
        if (c < 8 || c < b + 6 * Math.sqrt(b + 1)) continue;
        let id = -1;
        for (const q of [x > 0 ? k - 1 : -1, y > 0 ? k - nx : -1]) {
          if (q >= 0 && cur[q] >= 0) id = id < 0 ? this.find(cur[q]) : this.union(id, cur[q]);
        }
        if (prev[k] >= 0) id = id < 0 ? this.find(prev[k]) : this.union(id, prev[k]);
        if (id < 0) {
          id = this.parent.length;
          this.parent.push(id);
          this.st.push([0, 0, 0, 0, 0, 0, z, z, 0]);
        }
        const w = c - b, s = this.st[id];
        s[0] += w; s[1] += w * x; s[2] += w * y; s[3] += w * (z + 0.5); s[4] += 1;
        if (c > s[5]) s[5] = c;
        if (z < s[6]) s[6] = z;
        if (z > s[7]) s[7] = z;
        if (near[k]) s[8] = 1;
        cur[k] = id;
      }
    }
    this.prev.set(cur);
    this.z++;
  }

  peaks(minPixels = 4, minNet = 150) {
    const out = [];
    for (let i = 0; i < this.parent.length; i++) {
      if (this.find(i) !== i) continue;
      const [w, wx, wy, wz, n, max, z0, z1, touches] = this.st[i];
      if (n < minPixels || w < minNet) continue;
      out.push({ x: wx / w, y: wy / w, z: wz / w, net: w, npix: n, max, frames: z1 - z0 + 1, touchesMask: !!touches });
    }
    return out;
  }
}

// ---- Geometry refinement ------------------------------------------------------------

const MT = transpose(M_CRYSALIS);
const W_ANGLE = 0.1; // 0.1 deg of scan-angle residual weighted like 1 pixel

/** Fraction of peaks within `tol` of integer indices for a UB and instrument model. */
export function indexPeaks(G, ub, peaks, tol = 0.15) {
  const inv = mul(inv3(ub), M_CRYSALIS);
  return peaks.map((p) => {
    const R = gonio(G, p.angles[0], p.angles[2], p.angles[3]);
    const P = pixelToLab(G, p.x, p.y), c = samplePosition(G, R);
    const D = [P[0] - c[0], P[1] - c[1], P[2] - c[2]], n = Math.hypot(...D);
    const x = [D[0] / n, D[1] / n, D[2] / n + 1];
    const h = mulv(inv, mulv(transpose(R), x));
    const dev = Math.max(...h.map((v) => Math.abs(v - Math.round(v))));
    return { h, hkl: h.map(Math.round), ok: dev < tol, dev };
  });
}

/**
 * Refine the instrument model and orientation against indexed peaks (each with x, y,
 * run, angles = header [omega, theta, kappa, phi] at the peak, axis = scan axis).
 * Levels as in docs/METHOD.md: L3 (orientation, cell, beam centre, d1, scan zero), L6
 * (+ crystal offset, kappa zero), L10 (goniometer fixed, one orientation per run, drift
 * knots for runs that still disagree). The detector distance stays at the header's
 * calibrated value: the peaks fix only cell / distance, so refining it would move the
 * absolute cell scale without improving the fit.
 */
export function refineGeometry(peaks, g0, ub0, { system = 'triclinic', onProgress = () => {} } = {}) {
  const sys = CELL_SYSTEMS[system];
  const cell0 = cellFromUB(ub0, g0.wavelength);
  const U0 = polarU(ub0, cell0, g0.wavelength);
  const cellStart = sys.names.map((n) => cell0[['a', 'b', 'c', 'alpha', 'beta', 'gamma'].indexOf(n)]);
  if (system === 'hexagonal' || system === 'tetragonal') cellStart[0] = 0.5 * (cell0[0] + cell0[1]);
  const runs = [...new Set(peaks.map((p) => p.run))].sort((a, b) => a - b);
  const counts = new Map(runs.map((r) => [r, peaks.filter((p) => p.run === r).length]));
  const ref = runs.reduce((a, b) => (counts.get(b) > counts.get(a) ? b : a), runs[0]);

  const build = (spec, p) => {
    const g = { ...spec.g };
    let o = 0;
    for (const k of spec.free) g[k] = p[o++];
    const rx = p[o++], ry = p[o++], rz = p[o++];
    const cellv = sys.cell(Array.from(p.slice(o, o + sys.names.length)));
    o += sys.names.length;
    const ub = mul(mul(smallRot([rx, ry, rz]), U0), bMatrix(...cellv)).map((v) => v * g.wavelength);
    const perRun = new Map(), drift = new Map();
    for (const r of spec.perRun) { perRun.set(r, [p[o], p[o + 1], p[o + 2]]); o += 3; }
    for (const [r, knots] of spec.drift) {
      const vals = [];
      for (let k = 0; k < knots.length; k++) { vals.push([p[o], p[o + 1], p[o + 2]]); o += 3; }
      drift.set(r, { knots, vals });
    }
    return { g, ub, cell: cellv, perRun, drift };
  };
  const delta = (m, run, angle) => {
    if (m.drift.has(run)) {
      const { knots, vals } = m.drift.get(run);
      const out = [0, 0, 0];
      const t = Math.min(Math.max(angle, knots[0]), knots.at(-1));
      let k = 0;
      while (k < knots.length - 2 && t > knots[k + 1]) k++;
      const f = (t - knots[k]) / (knots[k + 1] - knots[k] || 1);
      for (let c = 0; c < 3; c++) out[c] = vals[k][c] + f * (vals[k + 1][c] - vals[k][c]);
      return out;
    }
    return m.perRun.get(run) ?? null;
  };
  const resid = (m, set) => {
    const G = prepare(m.g), base = mul(MT, m.ub);
    const out = new Float64Array(3 * set.length);
    set.forEach((pk, n) => {
      let xC = mulv(base, pk.hkl);
      const d = delta(m, pk.run, pk.angles[pk.axis]);
      if (d) xC = mulv(smallRot(d), xC);
      const pr = predict(G, xC, pk.angles, pk.axis, pk.angles[pk.axis]);
      if (!pr) return;
      out[3 * n] = pk.x - pr.i;
      out[3 * n + 1] = pk.y - pr.j;
      out[3 * n + 2] = (pk.angles[pk.axis] - pr.angle) / W_ANGLE;
    });
    return out;
  };
  const STEP = { ox: 0.01, oy: 0.01, distance: 0.001, d1: 0.001, omegaOffset: 0.001, phiOffset: 0.001, kappaOffset: 0.001, tx: 1e-3, ty: 1e-3, tz: 1e-3 };
  const fit = (spec, set, label) => {
    const p0 = [...spec.free.map((k) => spec.g[k]), 0, 0, 0, ...cellStart];
    for (const r of spec.perRun) p0.push(...(spec.init?.perRun?.get(r) ?? [0, 0, 0]));
    for (const [r, knots] of spec.drift) for (let k = 0; k < knots.length; k++) p0.push(...(spec.init?.perRun?.get(r) ?? [0, 0, 0]));
    const steps = [...spec.free.map((k) => STEP[k]), 1e-3, 1e-3, 1e-3, ...cellStart.map((v, i) => (sys.names[i].length > 1 ? 1e-3 : v * 1e-5))];
    while (steps.length < p0.length) steps.push(1e-3);
    onProgress(label);
    const res = levenbergMarquardt((p) => resid(build(spec, p), set), p0, { steps });
    const m = build(spec, res.p);
    return { m, res, stats: residualStats(resid(m, set), set) };
  };

  const report = { referenceRun: ref, peaksPerRun: Object.fromEntries(counts), system, levels: {} };
  const det = ['ox', 'oy']; // distance fixed (see above)
  const scan0 = peaks[0]?.axis === 3 ? 'phiOffset' : 'omegaOffset';
  const L3 = fit({ g: g0, free: [...det, 'd1', scan0], perRun: [], drift: new Map() }, peaks, 'Refining the detector and orientation (L3)');
  report.levels.L3 = L3.stats;
  const L6 = fit({ g: L3.m.g, free: [...det, 'd1', scan0, 'tx', 'ty', 'tz', 'kappaOffset'], perRun: [], drift: new Map() }, peaks,
    'Adding the crystal offset and kappa zero (L6)');
  report.levels.L6 = L6.stats;
  const free10 = [...det, 'd1', 'tx', 'ty', 'tz'];
  const others = runs.filter((r) => r !== ref);
  const L10a = fit({ g: L6.m.g, free: free10, perRun: others, drift: new Map() }, peaks, 'One orientation per run (L10)');
  const rmsAngle = runs.map((r) => L10a.stats.perRun[r].rmsAngle);
  const med = median(rmsAngle);
  const driftRuns = runs.filter((r, i) => rmsAngle[i] > Math.max(0.15, 2 * med));
  report.levels.L10a = { ...L10a.stats, driftRuns };
  let final = L10a;
  if (driftRuns.length) {
    const drift = new Map(driftRuns.map((r) => {
      const a = peaks.filter((p) => p.run === r).map((p) => p.angles[p.axis]);
      const lo = Math.floor(Math.min(...a)), hi = Math.ceil(Math.max(...a));
      return [r, Array.from({ length: 6 }, (_, k) => lo + (hi - lo) * k / 5)];
    }));
    final = fit({ g: L6.m.g, free: free10, perRun: others.filter((r) => !drift.has(r)), drift, init: { perRun: L10a.m.perRun } },
      peaks, 'Drift within runs (L10)');
  }
  report.levels.L10 = final.stats;
  return { model: final.m, report, cellStart: cell0 };
}

function residualStats(r, set) {
  const pick = (idx) => {
    const dx = [], dy = [], da = [];
    for (const n of idx) { dx.push(r[3 * n]); dy.push(r[3 * n + 1]); da.push(r[3 * n + 2] * W_ANGLE); }
    const rms = (a) => Math.sqrt(a.reduce((s, v) => s + v * v, 0) / Math.max(1, a.length));
    const med = (a) => median(a.map(Math.abs));
    return { n: idx.length, rmsX: rms(dx), rmsY: rms(dy), rmsAngle: rms(da), medX: med(dx), medY: med(dy), medAngle: med(da) };
  };
  const all = pick(set.map((_, n) => n));
  const perRun = {};
  for (const run of new Set(set.map((p) => p.run))) perRun[run] = pick(set.map((p, n) => (p.run === run ? n : -1)).filter((n) => n >= 0));
  return { ...all, perRun };
}

// ---- Gridding ---------------------------------------------------------------------------

/** Grid with voxel centres at origin + k*step covering [min, max] on each axis (H, K, L). */
export function makeGrid(min, max, step, origin = [0, 0, 0]) {
  const lo = min.map((v, i) => origin[i] + Math.floor((v - origin[i]) / step + 1e-9) * step);
  const shape = max.map((v, i) => Math.floor((v - lo[i]) / step + 1e-9) + 1);
  return { min: lo, step, shape }; // shape = [nH, nK, nL]
}

/**
 * Accumulators for a grid: counts, Poisson variance, normalization weight, pixel-frames, and
 * with `background` the air-scatter counts expected in the same pixel-frames (B).
 */
export function makeAccumulators(grid, { background = false } = {}) {
  const n = grid.shape[0] * grid.shape[1] * grid.shape[2];
  try {
    const acc = { S: new Float64Array(n), E2: new Float32Array(n), W: new Float64Array(n), N: new Uint32Array(n) };
    if (background) acc.B = new Float64Array(n);
    return acc;
  } catch {
    throw new Error(`Could not allocate ${(n * (background ? 32 : 24) / 1e9).toFixed(2)} GB for ${n.toLocaleString()} voxels; use a larger voxel or a smaller range.`);
  }
}

/**
 * Bin one frame: for every listed pixel p (lab scattering vector x[3p..3p+2], counts c[p],
 * normalization weight w[p]) and every sub-sample matrix A[k] (H = A x), add its share of
 * the counts to the voxel it falls in. Shares landing in the same voxel are summed
 * before squaring for the variance, so split counts are not treated as independent.
 */
export function accumulateFrame(acc, grid, x, c, w, A) {
  const vox = new Int32Array(c.length);
  applyFrame(acc, vox, c, w, indexFrame(grid, x, A, vox), A.length);
}

/**
 * The voxels of one frame's pixels: vox[p] is the voxel index when all sub-samples of pixel p
 * fall in one voxel (-1 when that is outside the grid), else -2, and the pixel's nsub voxel
 * indices (-1 outside) follow, in pixel order, in the returned array.
 */
export function indexFrame(grid, x, A, vox) {
  const [nH, nK, nL] = grid.shape, step = grid.step, [h0, k0, l0] = grid.min;
  const nsub = A.length, np = vox.length;
  const a = new Float64Array(9 * nsub); // all sub-sample matrices in one typed array
  for (let k = 0; k < nsub; k++) a.set(A[k], 9 * k);
  let split = new Int32Array(1 << 16), ns = 0;
  for (let p = 0; p < np; p++) {
    const x0 = x[3 * p], x1 = x[3 * p + 1], x2 = x[3 * p + 2];
    if (ns + nsub > split.length) { const grown = new Int32Array(2 * split.length); grown.set(split); split = grown; }
    let same = true, first = 0;
    for (let k = 0, o = 0; k < nsub; k++, o += 9) {
      const ih = Math.floor(((a[o] * x0 + a[o + 1] * x1 + a[o + 2] * x2) - h0) / step + 0.5);
      const ik = Math.floor(((a[o + 3] * x0 + a[o + 4] * x1 + a[o + 5] * x2) - k0) / step + 0.5);
      const il = Math.floor(((a[o + 6] * x0 + a[o + 7] * x1 + a[o + 8] * x2) - l0) / step + 0.5);
      const v = ih < 0 || ih >= nH || ik < 0 || ik >= nK || il < 0 || il >= nL ? -1 : (il * nK + ik) * nH + ih;
      if (k === 0) first = v;
      else if (v !== first) same = false;
      split[ns + k] = v;
    }
    if (same) vox[p] = first; // the usual case: the whole frame sweep stays in one voxel
    else { vox[p] = -2; ns += nsub; }
  }
  return split.slice(0, ns);
}

/**
 * A split pixel's voxels (its nsub entries from split[o]): fn(u, m) for each voxel u inside the
 * grid, at its first occurrence, with the number m of sub-samples that fell in it.
 */
function splitTerms(split, o, nsub, fn) {
  for (let k = 0; k < nsub; k++) {
    const u = split[o + k];
    if (u < 0) continue;
    let seen = false;
    for (let q = 0; q < k; q++) if (split[o + q] === u) { seen = true; break; }
    if (seen) continue;
    let m = 1;
    for (let q = k + 1; q < nsub; q++) if (split[o + q] === u) m++;
    fn(u, m);
  }
}

/** Add one frame's contributions (indexFrame) to the accumulators, in pixel order. */
export function applyFrame(acc, vox, c, w, split, nsub, b = null) {
  const { S, E2, W, N, B } = acc, np = vox.length;
  let o = 0;
  for (let p = 0; p < np; p++) {
    const v = vox[p];
    if (v >= 0) { const cp = c[p]; S[v] += cp; E2[v] += cp; W[v] += w[p]; N[v] += 1; if (b) B[v] += b[p]; continue; }
    if (v === -1) continue;
    const cp = c[p], wp = w[p], bp = b ? b[p] : 0;
    splitTerms(split, o, nsub, (u, m) => {
      const f = m / nsub;
      S[u] += f * cp;
      E2[u] += f * f * cp;
      W[u] += f * wp;
      N[u] += 1;
      if (b) B[u] += f * bp;
    });
    o += nsub;
  }
}

/**
 * One frame's contributions as records grouped by voxel block (for the parallel path): voxel
 * v[r] gets s[r] counts, e[r] variance, w[r] weight and (with `b`) b[r] background counts,
 * exactly the terms applyFrame adds.
 * Records are ordered by block of 4096 voxels (a stable counting sort), so each voxel still
 * receives its terms in pixel order, as in applyFrame: the sums are bit-identical, and adding
 * the records stays within a cache-sized part of the accumulators at a time.
 */
export function frameRecords(vox, c, w, split, nsub, nvox, b = null) {
  const np = vox.length, start = new Int32Array((nvox >>> 12) + 2);
  // count the terms per block ...
  for (let p = 0, o = 0; p < np; p++) {
    const v = vox[p];
    if (v >= 0) start[(v >>> 12) + 1]++;
    else if (v === -2) { splitTerms(split, o, nsub, (u) => { start[(u >>> 12) + 1]++; }); o += nsub; }
  }
  for (let b = 1; b < start.length; b++) start[b] += start[b - 1];
  // ... then place them in pixel order, as applyFrame adds them (stable within each block)
  const n = start[start.length - 1];
  const out = { v: new Int32Array(n), s: new Float64Array(n), e: new Float64Array(n), w: new Float64Array(n), b: b ? new Float64Array(n) : null };
  const put = (u, sv, ev, wv, bv) => {
    const j = start[u >>> 12]++;
    out.v[j] = u; out.s[j] = sv; out.e[j] = ev; out.w[j] = wv;
    if (b) out.b[j] = bv;
  };
  for (let p = 0, o = 0; p < np; p++) {
    const v = vox[p];
    if (v >= 0) put(v, c[p], c[p], w[p], b ? b[p] : 0);
    else if (v === -2) {
      const cp = c[p], wp = w[p], bp = b ? b[p] : 0;
      splitTerms(split, o, nsub, (u, m) => { const f = m / nsub; put(u, f * cp, f * f * cp, f * wp, f * bp); });
      o += nsub;
    }
  }
  return out;
}

/** Add records (frameRecords) to the accumulators. */
export function applyRecords(acc, rec) {
  const { S, E2, W, N, B } = acc, { v, s, e, w, b } = rec, n = v.length;
  for (let r = 0; r < n; r++) { const u = v[r]; S[u] += s[r]; E2[u] += e[r]; W[u] += w[r]; N[u] += 1; }
  if (b) for (let r = 0; r < n; r++) B[v[r]] += b[r];
}

// ---- Per-frame work of each pass, shared by the serial path and rigaku-map-worker.js ----

/** Pass 1, one frame: check and decode it, and add its counts to the per-pixel sums. */
export function sumFrame(fmt, buf, h, img, sums, run, label) {
  const { nx, ny, sum, counted } = sums;
  if (h.nx !== nx || h.ny !== ny) throw new Error(`frame ${label}: ${h.nx} x ${h.ny} pixels, expected ${nx} x ${ny}`);
  if (h.scanAxis !== 0 && h.scanAxis !== 3) throw new Error(`frame ${label}: scan axis ${h.scanAxis} is not omega or phi`);
  fmt.decodeTY6(buf, h, img);
  if (!fmt.checkStats(h, img)) throw new Error(`frame ${label}: decoded pixels do not match the header statistics`);
  let rs = sums.runSum.get(run);
  if (!rs) sums.runSum.set(run, (rs = new Float64Array(nx * ny)));
  for (let k = 0; k < img.length; k++) { const v = img[k]; sum[k] += v; rs[k] += v; if (v > 0) counted[k] = 1; }
}

export const makeSums = (nx, ny) => ({ nx, ny, sum: new Float64Array(nx * ny), counted: new Uint8Array(nx * ny), runSum: new Map() });

/** Merge per-worker pass-1 sums (integer counts, so the order of addition does not matter). */
export function mergeSums(into, part) {
  for (let k = 0; k < into.sum.length; k++) { into.sum[k] += part.sum[k]; into.counted[k] |= part.counted[k]; }
  for (const [run, rs] of part.runSum) {
    const t = into.runSum.get(run);
    if (!t) into.runSum.set(run, rs);
    else for (let k = 0; k < t.length; k++) t[k] += rs[k];
  }
}

/** Pass 2, one run: its Bragg peaks, from its frames in order. */
export function runPeaks(fmt, bufs, hs, mask, bg, nx, ny) {
  const ps = new PeakSearch(nx, ny, mask, bg), img = new Int32Array(nx * ny);
  for (let k = 0; k < bufs.length; k++) { fmt.decodeTY6(bufs[k], hs[k], img); ps.addFrame(img); }
  return ps.peaks();
}

/**
 * Pass 3 context from plain data that a worker can receive: setup = { nx, ny, unmasked
 * (Int32Array of pixel indices), model ({ g, ub, perRun, drift }), T, nSub, sap, mono }.
 */
export function mapContext(setup) {
  const { nx, ny, unmasked, model } = setup, np = unmasked.length;
  const G = prepare(model.g);
  const lab = new Float64Array(3 * np);
  for (let p = 0; p < np; p++) { const k = unmasked[p]; lab.set(pixelToLab(G, k % nx, Math.floor(k / nx)), 3 * p); }
  const px2 = model.g.pixelMM ** 2;
  const deltaOf = (run, angle) => {
    if (model.drift.has(run)) {
      const { knots, vals } = model.drift.get(run);
      const t = Math.min(Math.max(angle, knots[0]), knots.at(-1));
      let k = 0;
      while (k < knots.length - 2 && t > knots[k + 1]) k++;
      const f = (t - knots[k]) / (knots[k + 1] - knots[k] || 1);
      return [0, 1, 2].map((q) => vals[k][q] + f * (vals[k + 1][q] - vals[k][q]));
    }
    return model.perRun.get(run) ?? [0, 0, 0];
  };
  return {
    ...setup, G, np, lab, normal: G.normal, px2, omegaRef: px2 / model.g.distance ** 2, pol: polarizationModel(setup.mono),
    UBinvM: mul(inv3(model.ub), M_CRYSALIS), deltaOf, img: new Int32Array(nx * ny),
  };
}

/**
 * Pass 3, one frame: decode it, fill out.x (lab scattering vectors), out.c (counts) and out.w
 * (normalization weights) for the unmasked pixels, and return the sub-sample matrices A.
 */
export function mapFrame(ctx, fmt, buf, h, run, out) {
  const { G, np, lab, unmasked, img, normal, px2, omegaRef, pol, sap, nSub, T, UBinvM, deltaOf } = ctx;
  const { x, c, w } = out;
  fmt.decodeTY6(buf, h, img);
  const axis = h.scanAxis, a0 = h.start[axis], a1 = h.end[axis];
  const mid = h.start.slice(0, 4);
  mid[axis] = 0.5 * (a0 + a1);
  const cpos = samplePosition(G, gonio(G, mid[0], mid[2], mid[3]));
  for (let p = 0; p < np; p++) {
    const Dx = lab[3 * p] - cpos[0], Dy = lab[3 * p + 1] - cpos[1], Dz = lab[3 * p + 2] - cpos[2];
    const r = Math.sqrt(Dx * Dx + Dy * Dy + Dz * Dz), sx = Dx / r, sy = Dy / r, sz = Dz / r;
    x[3 * p] = sx; x[3 * p + 1] = sy; x[3 * p + 2] = sz + 1;
    c[p] = img[unmasked[p]];
    if (sap) {
      const dOmega = px2 * (sx * normal[0] + sy * normal[1] + sz * normal[2]) / (r * r);
      w[p] = h.exposure * dOmega * pol(sx, sy, sz) / omegaRef;
    } else w[p] = h.exposure;
  }
  const A = [];
  for (let k = 0; k < nSub; k++) {
    const ang = mid.slice();
    ang[axis] = a0 + (k + 0.5) / nSub * (a1 - a0);
    const R = gonio(G, ang[0], ang[2], ang[3]);
    A.push(mul(mul(mul(T, UBinvM), transpose(smallRot(deltaOf(run, ang[axis])))), transpose(R)));
  }
  return A;
}

// Graphite monochromator polarization (see docs/METHOD.md): sigma along lab X when the
// monochromator's scattering plane is CrysAlis E1E3 (lab Y-Z), along Y for E1E2.
export function polarizationFactor(s, mono) {
  return polarizationModel(mono)(s[0], s[1], s[2]);
}

/** polarizationFactor as a function of the unit vector's components, constants worked out once. */
export function polarizationModel(mono) {
  if (!mono || !(mono.theta > 0)) return (sx, sy, sz) => 0.5 * (1 + (sz * sz)); // unpolarized: (1 + cos^2 2theta) / 2
  const c2 = Math.cos(2 * mono.theta * DEG) ** 2;
  if (/E1E2/i.test(mono.plane ?? '')) return (sx, sy) => ((1 - sy * sy) + c2 * (1 - sx * sx)) / (1 + c2);
  return (sx, sy) => ((1 - sx * sx) + c2 * (1 - sy * sy)) / (1 + c2);
}

// ---- Mantid MDHistoWorkspace (SaveMD version 2) ------------------------------------------

/**
 * Write an MDHistoWorkspace as SaveMD v2 lays it out (checked against a file written by
 * Mantid 6.16 and reloaded with LoadMD): /MDHistoWorkspace/{data/{signal, errors_squared,
 * mask, num_events, D0, D1, D2}, experiment0/{instrument, logs, sample/oriented_lattice}}.
 * `signal`, `errors2`, `events` are Float64Arrays in (L, K, H) order (H fastest).
 */
export function writeMantidMD(h5wasm, path, { shape, edges, signal, errors2, events, ub, cell, logs = {}, title = '' }) {
  const f = new h5wasm.File(path, 'w');
  const str = (g, name, value) => {
    const v = String(value) || ' ';
    g.create_dataset({ name, data: [v], shape: [1], dtype: `A${new TextEncoder().encode(v).length}` });
  };
  try {
    f.create_attribute('NX_class', 'NXroot');
    const root = f.create_group('MDHistoWorkspace');
    root.create_attribute('NX_class', 'NXentry');
    root.create_attribute('QConvention', 'Crystallography');
    root.create_attribute('SaveMDVersion', 2, null, '<i4');
    root.create_dataset({ name: 'coordinate_system', data: Uint32Array.of(3) });
    root.create_dataset({ name: 'visual_normalization', data: Uint32Array.of(0) });
    if (title) str(root, 'title', title);
    const data = root.create_group('data');
    data.create_attribute('NX_class', 'NXdata');
    ['[H,0,0]', '[0,K,0]', '[0,0,L]'].forEach((long, i) => {
      const d = data.create_dataset({ name: `D${i}`, data: Float64Array.from(edges[i]) });
      d.create_attribute('frame', 'HKL');
      d.create_attribute('long_name', long);
      d.create_attribute('units', 'r.l.u.');
    });
    const chunks = [1, shape[1], shape[2]];
    const big = (name, arr) => data.create_dataset({ name, data: arr, shape, chunks, compression: 'gzip', compression_opts: 4 });
    const sig = big('signal', signal);
    sig.create_attribute('axes', 'D2:D1:D0');
    sig.create_attribute('signal', 1, null, '<i4');
    big('errors_squared', errors2);
    big('mask', new Int8Array(signal.length));
    big('num_events', events);
    const exp = root.create_group('experiment0');
    exp.create_attribute('NX_class', 'NXgroup');
    exp.create_attribute('version', 1, null, '<i4');
    const inst = exp.create_group('instrument');
    inst.create_attribute('NX_class', 'NXinstrument');
    inst.create_attribute('version', 1, null, '<i4');
    str(inst, 'name', ' ');
    const pmap = inst.create_group('instrument_parameter_map');
    pmap.create_attribute('NX_class', 'NXnote');
    pmap.create_attribute('version', 1, null, '<i4');
    for (const [k, v] of [['author', ' '], ['data', ' '], ['date', new Date().toISOString()], ['description', 'A string representation of the parameter map.'], ['type', 'text/plain']]) str(pmap, k, v);
    const xml = inst.create_group('instrument_xml');
    xml.create_attribute('NX_class', 'NXnote');
    for (const [k, v] of [['data', ' '], ['description', 'XML contents of the instrument IDF file.'], ['type', 'text/xml']]) str(xml, k, v);
    const lg = exp.create_group('logs');
    lg.create_attribute('NX_class', 'NXgroup');
    lg.create_attribute('version', 1, null, '<i4');
    const log = (name, value, units = ' ') => {
      const g = lg.create_group(name);
      g.create_attribute('NX_class', 'NXlog');
      const d = typeof value === 'string' ? (str(g, 'value', value), g.get('value'))
        : g.create_dataset({ name: 'value', data: Float64Array.from([].concat(value)) });
      d.create_attribute('units', units);
    };
    log('W_MATRIX', [1, 0, 0, 0, 1, 0, 0, 0, 1]);
    const gon = lg.create_group('goniometer');
    gon.create_attribute('NX_class', 'NXpositioner');
    gon.create_attribute('version', 1, null, '<i4');
    gon.create_dataset({ name: 'num_axes', data: Int32Array.of(0) });
    gon.create_dataset({ name: 'rotation_matrix', data: Float64Array.of(1, 0, 0, 0, 1, 0, 0, 0, 1) });
    for (const [name, v] of Object.entries(logs)) {
      if (v && typeof v === 'object' && 'value' in v) log(name, v.value, v.units ?? ' ');
      else log(name, typeof v === 'number' ? v : String(v));
    }
    const sample = exp.create_group('sample');
    sample.create_attribute('NX_class', 'NXsample');
    sample.create_attribute('name', ' ');
    sample.create_attribute('name_empty', 1, null, '<i4');
    sample.create_attribute('shape_xml', '<type name="userShape">  </type>');
    sample.create_attribute('version', 1, null, '<i4');
    for (const k of ['geom_height', 'geom_thickness', 'geom_width']) sample.create_dataset({ name: k, data: Float64Array.of(0) });
    sample.create_dataset({ name: 'geom_id', data: Int32Array.of(0) });
    sample.create_dataset({ name: 'num_oriented_lattice', data: Int32Array.of(1) });
    sample.create_dataset({ name: 'num_other_samples', data: Int32Array.of(0) });
    const mat = sample.create_group('material');
    mat.create_attribute('NX_class', 'NXdata');
    mat.create_attribute('formulaStyle', 'empty');
    mat.create_attribute('name', ' ');
    mat.create_attribute('version', 2, null, '<i4');
    for (const [k, v] of [['number_density', 0], ['packing_fraction', 1], ['pressure', 0], ['temperature', 0]]) mat.create_dataset({ name: k, data: Float64Array.of(v) });
    const ol = sample.create_group('oriented_lattice');
    ol.create_attribute('NX_class', 'NXcrystal');
    ol.create_dataset({ name: 'cross_term', data: Int32Array.of(0) });
    ol.create_dataset({ name: 'maximum_order', data: Int32Array.of(0) });
    for (const k of ['modulated_hkl_error', 'modulated_orientation_matrix']) ol.create_dataset({ name: k, data: new Float64Array(9), shape: [3, 3] });
    ol.create_dataset({ name: 'orientation_matrix', data: Float64Array.from(ub), shape: [3, 3] });
    ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].forEach((k, i) => {
      ol.create_dataset({ name: `unit_cell_${k}`, data: Float64Array.of(cell[i]) });
      ol.create_dataset({ name: `unit_cell_${k}_error`, data: Float64Array.of(0) });
    });
  } finally {
    f.close();
  }
}

// ---- The whole reduction --------------------------------------------------------------

const C_MANTID = [-1, 0, 0, 0, 1, 0, 0, 0, -1]; // lab (beam -Z) -> Mantid (beam +Z, Y up)

// Chip-boundary triplets of the HyPix-3000, found in the data sets of docs/METHOD.md (Rigaku reduction).
export const HYPIX3000_BOUNDARIES = {
  columns: [95, 96, 97, 192, 193, 194, 289, 290, 291, 386, 387, 388, 483, 484, 485, 580, 581, 582, 677, 678, 679],
  rows: [191, 192, 193],
};

/**
 * Reduce a CrysAlisPro experiment. `source` = { runs: Map(run -> [{ frame, read: async
 * () => ArrayBuffer }]), ubCandidates: { name: ub9 }, laue, monochromator, temperature,
 * label }, with parseRodHeader/decodeTY6 passed in `fmt`. Options:
 *   runs (array, default all), transform (row-major 3x3, rows = output basis vectors in
 *   units of the refined cell's: a'_i = sum_j T_ij a_j, so output indices are T x refined
 *   indices; determinant > 0), or multiplier (transform = multiplier x identity), or targetA
 *   (multiplier = round(targetA / a)); step (output r.l.u.), origin (a voxel centre), range ({min, max} in output r.l.u.,
 *   default: everything the detector reaches), nSub (sub-samples per frame, 5),
 *   normalization ('sap': solid angle + polarization, or 'rate'), refine (true),
 *   boundaries (extra chip-boundary lines).
 * `source.background` = { frames: [{ run, frame, read }], mode: 'static' | 'rotation', omegaBin }
 * adds a measured background as the accumulator B (finishVolume subtracts it): static, for air
 * scatter with the crystal and mount out of the beam; rotation, for an empty mount scanned like
 * the experiment (each sample run uses the background run with the same fixed angles).
 * `pool` (a WorkerPool of rigaku-map-worker.js) runs the per-frame work in parallel with
 * bit-identical results.
 * Returns { acc, grid, model, report } (see writeReduced to make the NeXus file).
 */
export async function reduceRigaku(source, fmt, opts = {}, onProgress = () => {}, pool = null) {
  const t0 = Date.now();
  const runIds = (opts.runs ?? [...source.runs.keys()]).filter((r) => source.runs.has(r));
  const frames = runIds.flatMap((r) => source.runs.get(r).map((f) => ({ ...f, run: r })));
  if (!frames.length) throw new Error('No frames to reduce.');
  const report = { label: source.label, runs: runIds, frames: frames.length, options: { ...opts, boundaries: undefined } };
  let done = 0;
  const bgFrames = source.background?.frames ?? [], bgMode = source.background?.mode ?? 'static';
  const bgWidth = source.background?.omegaBin ?? 5;
  if (bgFrames.length && !['static', 'rotation'].includes(bgMode)) throw new Error(`Unknown background mode "${bgMode}".`);
  if (bgFrames.length && !(bgWidth > 0)) throw new Error('The background bin width must be positive.');
  const total = 3 * frames.length + 2 * bgFrames.length;
  const tick = (label) => onProgress(label, ++done / total);

  report.threads = pool ? pool.size : 1;

  // Pass 1: sums, headers
  const first = fmt.parseRodHeader(await frames[0].read());
  const { nx, ny } = first;
  // Background frames: setup checked before any heavy work
  const bgIndex = bgFrames.length ? await backgroundIndex(fmt, bgFrames, first, () => tick('Checking the background frames')) : null;
  const sums = makeSums(nx, ny);
  const headers = new Array(frames.length);
  const img = new Int32Array(nx * ny);
  const label1 = 'Reading frames (1 of 3)';
  if (pool) {
    await pool.all({ type: 'sum-begin', nx, ny });
    await pool.ordered(frames.length, async (i) => {
      const f = frames[i], buf = await f.read();
      return { msg: { type: 'sum', buf, run: f.run, label: `${f.run}_${f.frame}` }, transfer: [buf] };
    }, (i, d) => { headers[i] = d.h; tick(label1); });
    for (const part of await pool.all({ type: 'sum-end' })) mergeSums(sums, part);
  } else {
    for (let i = 0; i < frames.length; i++) {
      const f = frames[i], buf = await f.read();
      const h = fmt.parseRodHeader(buf);
      sumFrame(fmt, buf, h, img, sums, f.run, `${f.run}_${f.frame}`);
      headers[i] = h;
      tick(label1);
    }
  }
  const { sum, counted, runSum } = sums;
  const maxCounts = Math.max(...headers.map((h) => h.stat.max));
  report.frameCheck = { decoded: frames.length, statsMatched: frames.length, maxPixelCounts: maxCounts,
    atOverflow: headers.filter((h) => h.stat.max >= h.overflowThreshold).length };
  // HyPix-3000 (775 x 385): 2 x 8 chips; their boundary lines are masked even where their
  // contrast is too weak to detect (it varies between data sets).
  const known = nx === 775 && ny === 385 && opts.knownBoundaries !== false ? HYPIX3000_BOUNDARIES : null;
  const boundaries = known || opts.boundaries
    ? { columns: [...(known?.columns ?? []), ...(opts.boundaries?.columns ?? [])], rows: [...(known?.rows ?? []), ...(opts.boundaries?.rows ?? [])] } : null;
  const { mask, info: maskInfo } = opts.mask ? { mask: opts.mask, info: { given: true } }
    : detectorMask(sum, counted, nx, ny, { boundaries });
  if (known) maskInfo.knownBoundaries = 'HyPix-3000';
  report.mask = maskInfo;

  // Pass 2: peaks (skipped when the geometry is given)
  const g0 = headerGeometry(first);
  if (opts.model) done += frames.length;
  const peaks = [];
  const peakRuns = opts.model ? [] : runIds, label2 = 'Finding Bragg peaks (2 of 3)';
  const runStart = []; // index of each run's first frame in frames / headers
  for (let i = 0, fi = 0; i < runIds.length; fi += source.runs.get(runIds[i]).length, i++) runStart.push(fi);
  const runFound = new Array(peakRuns.length);
  const background = (r) => Float64Array.from(runSum.get(r), (v) => v / source.runs.get(r).length);
  if (pool) {
    await pool.ordered(peakRuns.length, async (ri) => {
      const r = peakRuns[ri], list = source.runs.get(r), fi = runStart[runIds.indexOf(r)];
      const bufs = await Promise.all(list.map((f) => f.read())), bg = background(r);
      return { msg: { type: 'peaks', bufs, hs: headers.slice(fi, fi + list.length), mask, bg, nx, ny }, transfer: [...bufs, bg.buffer] };
    }, (ri, d) => {
      runFound[ri] = d.peaks;
      done += source.runs.get(peakRuns[ri]).length;
      onProgress(label2, done / total);
    }, 1);
  } else {
    for (let ri = 0; ri < peakRuns.length; ri++) {
      const r = peakRuns[ri], list = source.runs.get(r), fi = runStart[runIds.indexOf(r)];
      const ps = new PeakSearch(nx, ny, mask, background(r));
      for (let k = 0; k < list.length; k++) {
        fmt.decodeTY6(await list[k].read(), headers[fi + k], img);
        ps.addFrame(img);
        tick(label2);
      }
      runFound[ri] = ps.peaks();
    }
  }
  peakRuns.forEach((r, ri) => {
    const h0 = headers[runStart[runIds.indexOf(r)]], axis = h0.scanAxis, width = h0.end[axis] - h0.start[axis];
    for (const p of runFound[ri]) {
      const angles = h0.start.slice(0, 4);
      angles[axis] = h0.start[axis] + width * p.z;
      peaks.push({ ...p, run: r, axis, angles });
    }
  });
  const clean = peaks.filter((p) => !p.touchesMask);
  report.peaks = { found: peaks.length, clean: clean.length };

  // Geometry
  let model;
  if (opts.model) {
    model = opts.model;
    report.geometry = { refined: false, given: true };
  } else {
    const G0 = prepare(g0);
    const cands = Object.entries(source.ubCandidates ?? {}).filter(([, ub]) => ub && ub.length === 9);
    if (!cands.length) throw new Error('No orientation matrix found (crystal.ini, *.par or *_cracker.par).');
    const scores = cands.map(([name, ub]) => [name, indexPeaks(G0, ub, clean, 0.1).filter((x) => x.ok).length / Math.max(1, clean.length)]);
    scores.sort((a, b) => b[1] - a[1]);
    report.ubCandidates = Object.fromEntries(scores);
    const [startName, startScore] = scores[0];
    if (startScore < 0.5) throw new Error(`The orientation matrices index only ${(100 * startScore).toFixed(0)} % of the peaks; the geometry conventions may not fit this instrument.`);
    const ub0 = cands.find(([n]) => n === startName)[1];
    const indexed = indexPeaks(G0, ub0, clean, 0.15);
    const set = clean.map((p, i) => ({ ...p, hkl: indexed[i].hkl })).filter((_, i) => indexed[i].ok);
    report.startUB = startName;
    report.indexed = set.length;
    const startCell = cellFromUB(ub0, g0.wavelength);
    const system = crystalSystem(source.laue, startCell);
    if (opts.refine === false) {
      model = { g: g0, ub: ub0, cell: startCell, perRun: new Map(), drift: new Map() };
      report.geometry = { refined: false, system };
    } else {
      const t = Date.now();
      const out = refineGeometry(set, g0, ub0, { system, onProgress: (label) => onProgress(label, done / total) });
      model = out.model;
      report.geometry = { refined: true, ...out.report, seconds: (Date.now() - t) / 1000 };
    }
  }
  report.cell = model.cell;
  report.model = {
    geometry: model.g, ub: model.ub, cell: model.cell, perRun: Object.fromEntries(model.perRun),
    drift: Object.fromEntries([...model.drift].map(([r, d]) => [r, d])),
  };

  // Grid
  const G = prepare(model.g);
  // output indices = T x indices of the refined cell (T = m x identity for a multiple)
  const m = opts.multiplier ?? (opts.targetA ? Math.max(1, Math.round(opts.targetA / model.cell[0])) : 1);
  const T = opts.transform ? [...opts.transform] : [m, 0, 0, 0, m, 0, 0, 0, m];
  const detT = det3(T);
  if (!(detT > 1e-9)) throw new Error('The output-cell transformation must have a positive determinant.');
  const cellOut = transformCell(model.cell, T);
  const isMultiple = T.every((v, i) => (i % 4 === 0 ? v === T[0] : v === 0));
  const step = opts.step ?? 0.05 * Math.cbrt(detT);
  let maxX = 0;
  for (let k = 0; k < nx * ny; k++) {
    if (mask[k]) continue;
    const P = pixelToLab(G, k % nx, Math.floor(k / nx)), n = Math.hypot(...P);
    maxX = Math.max(maxX, Math.hypot(P[0] / n, P[1] / n, P[2] / n + 1));
  }
  const qmax = maxX / model.g.wavelength;
  const lengths = cellOut.slice(0, 3);
  const reach = lengths.map((a) => Math.ceil(qmax * a / step) * step);
  const range = opts.range ?? { min: reach.map((v) => -v), max: reach };
  const grid = makeGrid(range.min, range.max, step, opts.origin ?? [0, 0, 0]);
  const nvox = grid.shape[0] * grid.shape[1] * grid.shape[2];
  if (nvox > (opts.maxVoxels ?? 6e7)) throw new Error(`${nvox.toLocaleString()} voxels: use a larger voxel or a smaller range.`);
  report.grid = { ...grid, voxels: nvox, qmax, transform: T, multiplier: isMultiple ? T[0] : null, cell: cellOut };
  // Pass 3: gridding
  const unmasked = [];
  for (let k = 0; k < nx * ny; k++) if (!mask[k]) unmasked.push(k);

  // Background: the counts it puts in each sample frame's unmasked pixels (rate x exposure)
  let frameBackground = null;
  if (bgIndex) {
    const bgLabel = 'Reading the background frames', stats = { matched: {} };
    const relErrors = (counts) => { const r = []; for (const v of counts) if (v > 0) r.push(1 / Math.sqrt(v)); r.sort((a, b) => a - b); return r.length ? r[r.length >> 1] : null; };
    if (bgMode === 'static') {
      const model = await staticBackground(fmt, bgIndex, first, unmasked, () => tick(bgLabel));
      stats.medianRelativeError = relErrors(model.counts);
      stats.meanRate = model.rate.reduce((a, v) => a + v, 0) / model.rate.length;
      frameBackground = async (i) => frameBackgroundCounts(model, headers[i]);
    } else {
      // every sample run needs a background run with the same scan axis and fixed angles
      const missing = runIds.filter((r) => { const h = headers[runStart[runIds.indexOf(r)]]; return !matchBackgroundRun(bgIndex, h.scanAxis, h.start); });
      if (missing.length) {
        const desc = missing.map((r) => { const h = headers[runStart[runIds.indexOf(r)]]; return `run ${r} (${[0, 2, 3].filter((k) => k !== h.scanAxis).map((k) => `${ANGLES[k]} ${Number(h.start[k].toFixed(2))}`).join(', ')})`; });
        throw new Error(`No background run was measured like ${desc.join('; ')}.`);
      }
      let cached = null, rels = [], rateSum = 0, rateN = 0;
      frameBackground = async (i) => {
        const h = headers[i], bgRun = matchBackgroundRun(bgIndex, h.scanAxis, h.start);
        if (cached?.run !== bgRun.run) {
          cached = await rotationBackground(fmt, bgRun, first, unmasked, bgWidth, () => tick(bgLabel));
          stats.matched[frames[i].run] = bgRun.run;
          for (const c of cached.counts) { rels.push(relErrors(c)); }
          const tot = new Float64Array(unmasked.length);
          cached.counts.forEach((c) => { for (let p = 0; p < tot.length; p++) tot[p] += c[p]; });
          const ex = cached.exposure.reduce((a, v) => a + v, 0);
          for (let p = 0; p < tot.length; p++) { rateSum += tot[p] / ex; rateN++; }
          stats.binExposure = [Math.min(stats.binExposure?.[0] ?? Infinity, ...cached.exposure), Math.max(stats.binExposure?.[1] ?? 0, ...cached.exposure)];
          rels = rels.filter((v) => v !== null).sort((a, b) => a - b);
          stats.medianRelativeError = rels.length ? rels[rels.length >> 1] : null;
          stats.meanRate = rateSum / rateN;
        }
        stats.matched[frames[i].run] = bgRun.run;
        return frameBackgroundCounts(cached, h);
      };
    }
    report.background = {
      mode: bgMode, frames: bgIndex.frames, runs: bgIndex.runs.length, exposure: bgIndex.exposure, stats,
      ...(bgMode === 'rotation' ? { binWidth: bgWidth } : {}),
      method: bgMode === 'static'
        ? 'static: per-pixel rate = summed counts / summed exposure over all background frames; x each sample frame\'s exposure, accumulated as B; signal = (S - scale x B) / W'
        : `rotation-resolved: each sample run uses the background run with the same scan axis and fixed angles; per-pixel rate = summed counts / summed exposure in ${bgWidth} deg bins of the scan angle, interpolated linearly to each frame's mid angle; x the frame's exposure, accumulated as B; signal = (S - scale x B) / W`,
      errors: 'the errors are the Poisson errors of the measured counts; the background\'s own statistical error (median relative error per pixel and bin in stats) is not propagated',
    };
  }
  const acc = makeAccumulators(grid, { background: !!bgIndex });
  const nSub = opts.nSub ?? 5, label3 = 'Mapping pixels to HKL (3 of 3)';
  const setup = {
    nx, ny, unmasked: Int32Array.from(unmasked), nSub, T, sap: (opts.normalization ?? 'sap') === 'sap', mono: source.monochromator,
    model: { g: model.g, ub: model.ub, perRun: model.perRun, drift: model.drift },
  };
  const ctx = mapContext(setup), np = ctx.np;
  if (pool) {
    await pool.all({ type: 'map-init', setup });
    await pool.ordered(frames.length, async (i) => {
      const f = frames[i], buf = await f.read(), b = frameBackground ? await frameBackground(i) : null;
      return { msg: { type: 'map', buf, h: headers[i], run: f.run, grid, b }, transfer: [buf, ...(b ? [b.buffer] : [])] };
    }, (i, d) => { applyRecords(acc, d); tick(label3); });
  } else {
    const out = { x: new Float64Array(3 * np), c: new Int32Array(np), w: new Float64Array(np) };
    const vox = new Int32Array(np);
    for (let i = 0; i < frames.length; i++) {
      const A = mapFrame(ctx, fmt, await frames[i].read(), headers[i], frames[i].run, out);
      const b = frameBackground ? await frameBackground(i) : null;
      applyFrame(acc, vox, out.c, out.w, indexFrame(grid, out.x, A, vox), nSub, b);
      tick(label3);
    }
  }
  report.seconds = (Date.now() - t0) / 1000;
  report.normalization = setup.sap ? 'counts / (s x pixel solid angle x polarization) x Omega_ref' : 'counts / (s x pixel)';
  report.omegaRef = ctx.omegaRef;
  const ubMantid = mul(mul(C_MANTID, mul(MT, model.ub)), inv3(T)).map((v) => v / model.g.wavelength);
  return { acc, grid, model, report, ubMantid };
}

/**
 * Signal, errors^2 and pixel-frame counts in (L, K, H) order; NaN where nothing was
 * measured. With a background accumulator, signal = (S - backgroundScale x B) / W; the errors
 * stay the Poisson errors of the measured counts. `zeros` counts voxels with no measured
 * counts, `negative` (with a background) covered voxels whose signal is below zero. With
 * `inPlace`, the signal reuses acc.S and the counts acc.W (saving memory).
 */
export function finishVolume(acc, { inPlace = false, backgroundScale = 1 } = {}) {
  const n = acc.S.length, B = acc.B ?? null;
  const signal = inPlace ? acc.S : new Float64Array(n), errors2 = new Float64Array(n);
  let covered = 0, zeros = 0, negative = 0;
  for (let v = 0; v < n; v++) {
    const W = acc.W[v];
    if (W > 0) {
      const S = acc.S[v];
      if (S === 0) zeros++;
      signal[v] = (B ? S - backgroundScale * B[v] : S) / W;
      if (signal[v] < 0) negative++;
      errors2[v] = acc.E2[v] / (W * W);
      covered++;
    } else {
      signal[v] = NaN;
      errors2[v] = NaN;
    }
  }
  const events = inPlace ? acc.W : new Float64Array(n);
  for (let v = 0; v < n; v++) events[v] = errors2[v] === errors2[v] ? acc.N[v] : 0;
  return { signal, errors2, events, covered, zeros, ...(B ? { negative } : {}) };
}

// ---- Background: air scatter (static) or an empty mount scanned like the experiment ----------

/** How a frame's detector setup differs from the experiment's (`ref`, its first header). */
function setupDifferences(h, ref) {
  const near = (a, b, tol) => Math.abs(a - b) <= tol;
  return [
    (h.nx !== ref.nx || h.ny !== ref.ny) && `detector ${h.nx} x ${h.ny} pixels`,
    !near(h.distance, ref.distance, 1e-3) && `distance ${h.distance} mm, not ${ref.distance}`,
    !near(h.start[1], ref.start[1], 1e-3) && `2theta arm ${h.start[1]} deg, not ${ref.start[1]}`,
    !(h.pixelMM[0] === ref.pixelMM[0] && h.binning[0] === ref.binning[0] && h.binning[1] === ref.binning[1]) && 'different pixel size or binning',
    !near(h.wavelengths.alpha12, ref.wavelengths.alpha12, 1e-6) && `wavelength ${h.wavelengths.alpha12} A`,
  ].filter(Boolean);
}

const ANGLES = ['omega', 'theta', 'kappa', 'phi'];
const angleDiff = (a, b) => Math.abs(((a - b) % 360 + 540) % 360 - 180);

/**
 * Background frames, read once for their headers: every frame must have the experiment's
 * detector setup. Grouped by run, each with its scan axis, fixed angles and frames in scan order.
 */
export async function backgroundIndex(fmt, frames, ref, onFrame = () => {}) {
  const runs = new Map();
  let exposure = 0;
  for (let i = 0; i < frames.length; i++) {
    const f = frames[i], h = fmt.parseRodHeader(await f.read()), label = `background frame ${f.run ?? ''}_${f.frame ?? i + 1}`;
    const why = setupDifferences(h, ref);
    if (why.length) throw new Error(`${label}: not the experiment's setup (${why.join('; ')}).`);
    exposure += h.exposure;
    const r = f.run ?? 1;
    if (!runs.has(r)) runs.set(r, { run: r, axis: h.scanAxis, fixed: h.start.slice(0, 4), frames: [] });
    runs.get(r).frames.push({ f, h });
    onFrame();
  }
  for (const r of runs.values()) r.frames.sort((a, b) => a.h.start[r.axis] - b.h.start[r.axis]);
  if (!(exposure > 0)) throw new Error('The background frames have no exposure.');
  return { runs: [...runs.values()], exposure, frames: frames.length };
}

/**
 * The background run measured like a sample run: the same scan axis and the same fixed angles
 * (the goniometer angles that do not scan), within 0.01 deg. Null when there is none.
 */
export function matchBackgroundRun(index, axis, fixed) {
  return index.runs.find((r) => r.axis === axis && [0, 2, 3].every((k) => k === axis || angleDiff(r.fixed[k], fixed[k]) <= 0.01)) ?? null;
}

/** Decode and check one background frame; add its counts (unmasked pixels) to `into`. */
async function addBackgroundFrame(fmt, { f, h }, img, unmasked, into, label) {
  const buf = await f.read();
  fmt.decodeTY6(buf, h, img);
  if (!fmt.checkStats(h, img)) throw new Error(`${label}: decoded pixels do not match the header statistics`);
  for (let p = 0; p < unmasked.length; p++) into[p] += img[unmasked[p]];
}

/**
 * Static background (air scatter with the crystal and its mount out of the beam): for every
 * unmasked pixel, rate = summed counts / summed exposure over all the frames.
 */
export async function staticBackground(fmt, index, ref, unmasked, onFrame = () => {}) {
  const counts = new Float64Array(unmasked.length), img = new Int32Array(ref.nx * ref.ny);
  for (const r of index.runs) for (const fr of r.frames) { await addBackgroundFrame(fmt, fr, img, unmasked, counts, `background frame ${r.run}`); onFrame(); }
  return { mode: 'static', rate: Float64Array.from(counts, (v) => v / index.exposure), counts, exposure: index.exposure };
}

/**
 * Rotation-resolved background from one background run (an empty mount scanned like the
 * sample run): its frames are grouped into bins of `width` degrees of the scan angle, and each
 * bin gets rate = summed counts / summed exposure per unmasked pixel, at the exposure-weighted
 * mean angle of its frames (as mdx2 bins its background image series).
 */
export async function rotationBackground(fmt, run, ref, unmasked, width, onFrame = () => {}) {
  const img = new Int32Array(ref.nx * ref.ny), axis = run.axis, a0 = run.frames[0].h.start[axis];
  const bins = new Map();
  for (const fr of run.frames) {
    const mid = 0.5 * (fr.h.start[axis] + fr.h.end[axis]), j = Math.floor((mid - a0) / width + 1e-9);
    if (!bins.has(j)) bins.set(j, { counts: new Float64Array(unmasked.length), exposure: 0, angleSum: 0 });
    const b = bins.get(j);
    await addBackgroundFrame(fmt, fr, img, unmasked, b.counts, `background frame ${run.run}`);
    b.exposure += fr.h.exposure;
    b.angleSum += mid * fr.h.exposure;
    onFrame();
  }
  const list = [...bins.keys()].sort((x, y) => x - y).map((j) => bins.get(j));
  const first = run.frames[0].h, last = run.frames.at(-1).h;
  return {
    mode: 'rotation', axis, run: run.run, width, range: [first.start[axis], last.end[axis]],
    centres: list.map((b) => b.angleSum / b.exposure), exposure: list.map((b) => b.exposure),
    rates: list.map((b) => Float64Array.from(b.counts, (v) => v / b.exposure)), counts: list.map((b) => b.counts),
  };
}

/**
 * Background counts expected in one sample frame, per unmasked pixel: rate x the frame's
 * exposure. For a rotation-resolved background the rate is interpolated linearly in the scan
 * angle between bin centres (held at the end bins), at the frame's mid angle.
 */
export function frameBackgroundCounts(model, h) {
  if (model.mode === 'static') return Float64Array.from(model.rate, (v) => v * h.exposure);
  const axis = model.axis, mid = 0.5 * (h.start[axis] + h.end[axis]), c = model.centres, n = c.length;
  if (mid < model.range[0] - 1e-6 || mid > model.range[1] + 1e-6) {
    throw new Error(`The background run ${model.run} covers ${ANGLES[axis]} ${model.range[0]} to ${model.range[1]} deg, not ${mid}.`);
  }
  let r0 = model.rates[0], r1 = r0, u = 0;
  if (mid >= c[n - 1]) { r0 = r1 = model.rates[n - 1]; } else if (mid > c[0]) {
    let j = 0;
    while (c[j + 1] <= mid) j++;
    r0 = model.rates[j]; r1 = model.rates[j + 1]; u = (mid - c[j]) / (c[j + 1] - c[j]);
  }
  const out = new Float64Array(r0.length), t = h.exposure;
  for (let p = 0; p < out.length; p++) out[p] = ((1 - u) * r0[p] + u * r1[p]) * t;
  return out;
}

/** Bin edges of a grid, per axis (H, K, L). */
export const gridEdges = (grid) => grid.shape.map((n, i) => Float64Array.from({ length: n + 1 }, (_, k) => grid.min[i] + (k - 0.5) * grid.step));
