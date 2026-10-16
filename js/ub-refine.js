// Refine the UB matrix of an HKL volume against its own Bragg peaks (pure functions): find
// the peaks, index them with the volume's UB, and fit an orientation and a cell, which may
// be constrained to a crystal system, to their positions in Q.
//
// The volume's coordinates were computed with its UB (Mantid's convention, no 2π), so a
// peak found at indices h_obs was measured at Q = UB·h_obs, whatever error UB has. The
// refined UB′ = R·U₀·B(cell) minimizes Σ |UB′·h − Q|² over the peaks indexed h (integers),
// with U₀ the orientation of UB, R a small rotation and B the Busing–Levy matrix of the
// cell. Reindexing with UB′ (reindex.js) then puts the peaks on integer indices.

import { basisMatrix, uniformAxes } from './reindex.js';
import { CELL_SYSTEMS, DEG, bMatrix, cellFromUB, det3, inv3, levenbergMarquardt, mul, mulv, polarU, smallRot } from './rigaku-geometry.js';

/**
 * The crystal systems a cell can be constrained to, in standard settings (c unique for
 * tetragonal and hexagonal), with the parameters refined for each and how they expand
 * to a b c α β γ.
 */
export const UB_SYSTEMS = {
  cubic: { ...CELL_SYSTEMS.cubic, label: 'Cubic', short: 'cubic', free: 'a' },
  hexagonal: { ...CELL_SYSTEMS.hexagonal, label: 'Hexagonal, or trigonal on hexagonal axes', short: 'hexagonal', free: 'a, c' },
  rhombohedral: { names: ['a', 'alpha'], cell: ([a, al]) => [a, a, a, al, al, al], label: 'Trigonal on rhombohedral axes', short: 'rhombohedral', free: 'a, α' },
  tetragonal: { ...CELL_SYSTEMS.tetragonal, label: 'Tetragonal', short: 'tetragonal', free: 'a, c' },
  orthorhombic: { ...CELL_SYSTEMS.orthorhombic, label: 'Orthorhombic', short: 'orthorhombic', free: 'a, b, c' },
  monoclinic: { ...CELL_SYSTEMS.monoclinic, label: 'Monoclinic, b unique', short: 'monoclinic, b unique', free: 'a, b, c, β' },
  'monoclinic-c': { names: ['a', 'b', 'c', 'gamma'], cell: ([a, b, c, ga]) => [a, b, c, 90, 90, ga], label: 'Monoclinic, c unique', short: 'monoclinic, c unique', free: 'a, b, c, γ' },
  triclinic: { ...CELL_SYSTEMS.triclinic, label: 'Triclinic (unconstrained)', short: 'triclinic', free: 'a, b, c, α, β, γ' },
};

const NAMES = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'];
const mean = (...v) => v.reduce((s, x) => s + x, 0) / v.length;

/** The refined parameters of `system` that start from `cell` (equal lengths and angles averaged). */
export function startParams(system, cell) {
  const [a, b, c, al, be, ga] = cell;
  if (system === 'cubic') return [mean(a, b, c)];
  if (system === 'hexagonal' || system === 'tetragonal') return [mean(a, b), c];
  if (system === 'rhombohedral') return [mean(a, b, c), mean(al, be, ga)];
  return UB_SYSTEMS[system].names.map((n) => cell[NAMES.indexOf(n)]);
}

/** Whether `cell` has the metric of `system`: equal lengths within `length` (relative), angles within `angle` degrees. */
export function fitsSystem(system, cell, { length = 0.01, angle = 1 } = {}) {
  const ideal = UB_SYSTEMS[system].cell(startParams(system, cell));
  return ideal.every((x, i) => (i < 3 ? Math.abs(x - cell[i]) <= length * cell[i] : Math.abs(x - cell[i]) <= angle));
}

// Laue classes as the symmetry presets name them (trigonal ones on hexagonal axes).
const LAUE_SYSTEMS = [[/^m-3/, 'cubic'], [/^(6|-3)/, 'hexagonal'], [/^4/, 'tetragonal'], [/^mmm$/, 'orthorhombic'],
  [/^2\/m \(c/, 'monoclinic-c'], [/^2\/m/, 'monoclinic'], [/^-1$/, 'triclinic']];
// A cell must be this close to a system's metric to be refined in it: closer than a change of setting.
const LOOSE = { length: 0.05, angle: 5 };
const BY_SYMMETRY = ['cubic', 'hexagonal', 'rhombohedral', 'tetragonal', 'orthorhombic', 'monoclinic', 'monoclinic-c', 'triclinic'];

/**
 * The system to suggest for `cell`: that of the Laue class `laue` (a symmetry preset's
 * name) when the cell fits it, else the most symmetric one whose metric the cell has
 * within 1% and 1°. Returns { system, from: 'laue' | 'metric' }.
 */
export function suggestSystem(cell, laue = null) {
  const fromLaue = laue && LAUE_SYSTEMS.find(([re]) => re.test(laue))?.[1];
  if (fromLaue && fitsSystem(fromLaue, cell, LOOSE)) return { system: fromLaue, from: 'laue' };
  return { system: BY_SYMMETRY.find((s) => fitsSystem(s, cell)), from: 'metric' };
}

const median = (a) => {
  const s = Float64Array.from(a).sort();
  return s.length ? (s.length % 2 ? s[s.length >> 1] : 0.5 * (s[s.length / 2 - 1] + s[s.length / 2])) : NaN;
};

// Overlap of voxel i (spanning i ± 0.5) with the window c ± h.
const overlap = (i, c, h) => Math.max(0, Math.min(i + 0.5, c + h) - Math.max(i - 0.5, c - h));

/**
 * The centre of a peak near voxel `j` (display-axis voxel indices), in fractional voxel
 * indices: the intensity-weighted mean over a window of ± (`radius` - 1) voxels centred on
 * the estimate itself, each voxel weighted by its overlap with the window, less the
 * background, negative weights dropped, repeated until it moves less than 10⁻⁴ voxel. The
 * background is a plane through the box of ± `radius` voxels around the window's voxel:
 * the median of its surface, sloped along each axis by the difference of the medians of
 * the two opposite faces; a level background would leave the slope of diffuse scattering
 * or air scatter under the peak and pull its centre towards low |Q|. A window that follows
 * the estimate is symmetric about a symmetric peak's centre, so, unlike a box fixed on
 * voxels, it does not pull the centre towards the voxel grid. Null when the box leaves the
 * grid or holds an unmeasured voxel: next to a coverage gap, the bright voxels that detector
 * edges leave would bias the centre.
 */
function centroid(volume, n, j, radius) {
  const h = radius - 1;
  let c = j.slice(), out = null;
  for (let it = 0; it < 30; it++) {
    const centre = c.map(Math.round);
    if (centre.some((x, d) => x - radius < 0 || x + radius >= n[d])) return null;
    const shell = [], faces = [[], [], [], [], [], []]; // -x, +x, -y, +y, -z, +z
    for (let z = -radius; z <= radius; z++) {
      for (let y = -radius; y <= radius; y++) {
        for (let x = -radius; x <= radius; x++) {
          const v = volume[centre[0] + x + n[0] * (centre[1] + y + n[1] * (centre[2] + z))];
          if (v !== v) return null;
          if (Math.max(Math.abs(x), Math.abs(y), Math.abs(z)) < radius) continue;
          shell.push(v);
          if (Math.abs(x) === radius) faces[x > 0 ? 1 : 0].push(v);
          if (Math.abs(y) === radius) faces[y > 0 ? 3 : 2].push(v);
          if (Math.abs(z) === radius) faces[z > 0 ? 5 : 4].push(v);
        }
      }
    }
    const bg = median(shell), slope = [0, 1, 2].map((d) => (median(faces[2 * d + 1]) - median(faces[2 * d])) / (2 * radius));
    let sum = 0, sx = 0, sy = 0, sz = 0;
    for (let z = centre[2] - h; z <= centre[2] + h; z++) {
      const fz = overlap(z, c[2], h), bz = bg + slope[2] * (z - centre[2]);
      for (let y = centre[1] - h; y <= centre[1] + h; y++) {
        const fyz = fz * overlap(y, c[1], h), byz = bz + slope[1] * (y - centre[1]);
        for (let x = centre[0] - h; x <= centre[0] + h; x++) {
          const w = (volume[x + n[0] * (y + n[1] * z)] - byz - slope[0] * (x - centre[0])) * fyz * overlap(x, c[0], h);
          if (w > 0) { sum += w; sx += w * x; sy += w * y; sz += w * z; }
        }
      }
    }
    if (!(sum > 0)) return null;
    const next = [sx / sum, sy / sum, sz / sum];
    const moved = Math.max(...next.map((x, d) => Math.abs(x - c[d])));
    c = next;
    out = { c, intensity: sum, background: bg };
    if (moved < 1e-4) break;
  }
  return out;
}

/**
 * The half width at half maximum (voxels) of the strongest peaks along each display axis:
 * the median, over up to `count` local maxima `tops` (volume indices, strongest first) at
 * least 10 voxels apart and with every voxel within 2 measured (not the bright voxels at a
 * coverage edge), of the mean distance on either side at which the value falls to halfway
 * between the maximum and `level`, interpolated between voxels. NaN for an axis with no
 * such crossing.
 */
function halfWidths(volume, n, tops, level, count = 25) {
  const stride = [1, n[0], n[0] * n[1]], picked = [], widths = [[], [], []];
  const measured = (j) => {
    if (j.some((x, d) => x < 2 || x + 2 >= n[d])) return false;
    for (let z = -2; z <= 2; z++) {
      for (let y = -2; y <= 2; y++) {
        for (let x = -2; x <= 2; x++) if (Number.isNaN(volume[j[0] + x + n[0] * (j[1] + y + n[1] * (j[2] + z))])) return false;
      }
    }
    return true;
  };
  for (const v of tops) {
    if (picked.length >= count) break;
    const j = [v % n[0], Math.floor(v / n[0]) % n[1], Math.floor(v / stride[2])];
    if (picked.some((t) => Math.max(...t.map((x, d) => Math.abs(x - j[d]))) < 10) || !measured(j)) continue;
    picked.push(j);
    const half = 0.5 * (volume[v] + level);
    for (let d = 0; d < 3; d++) {
      const sides = [-1, 1].map((s) => {
        let prev = volume[v];
        for (let k = 1; j[d] + s * k >= 0 && j[d] + s * k < n[d]; k++) {
          const x = volume[v + s * k * stride[d]];
          if (x !== x) return NaN;
          if (x <= half) return k - 1 + (prev - half) / (prev - x);
          prev = x;
        }
        return NaN;
      });
      if (sides.every(Number.isFinite)) widths[d].push(0.5 * (sides[0] + sides[1]));
    }
  }
  return widths.map((w) => (w.length ? median(w) : NaN));
}

/**
 * Bragg peaks in `volume` (float32 in the order of `shape`, display axis d on volume axis
 * 2 - d, NaN where unmeasured), whose `dims` are HKL axes with uniform bins: the strongest
 * `maxPeaks` local maxima (over their 26 neighbours) above median + `threshold` × the
 * robust spread (1.4826 × MAD) of the measured voxels, each centred by `centroid` with a
 * box of ± `radius` voxels, and kept unless a stronger peak's centre lies within 2 ×
 * `radius` - 1 voxels of its maximum or its centre. Without `radius`, the box is
 * sized to the peaks: 1 + 2 × their largest half width along an axis, rounded up, from 3 to
 * 12 voxels. Returns { peaks: [{ h (indices), u (display coordinates), max, intensity,
 * background }], stats: { level, spread, threshold, maxima, halfWidth, radius, edge, near } }.
 */
export function findPeaks(volume, shape, dims, { threshold = 10, maxPeaks = 500, radius = null, onProgress = () => {} } = {}) {
  const axes = uniformAxes(dims), W = basisMatrix(dims);
  const [n0, n1, n2] = shape, n = [n2, n1, n0], plane = n1 * n2;
  // Level and spread of the measured voxels, from about a million of them.
  const stride = Math.max(1, Math.floor(volume.length / 1e6)), sample = [];
  for (let v = 0; v < volume.length; v += stride) if (volume[v] === volume[v]) sample.push(volume[v]);
  if (!sample.length) throw new Error('The volume has no measured voxels.');
  const level0 = median(sample);
  let spread = 1.4826 * median(sample.map((x) => Math.abs(x - level0)));
  if (!(spread > 0)) spread = Math.sqrt(sample.reduce((s, x) => s + (x - level0) ** 2, 0) / sample.length); // mostly equal values
  const level = level0 + threshold * spread;

  const offsets = [];
  for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) if (a || b || c) offsets.push(a * plane + b * n2 + c);
  const index = [], value = [];
  for (let i0 = 1; i0 < n0 - 1; i0++) {
    for (let i1 = 1; i1 < n1 - 1; i1++) {
      for (let i2 = 1, v = i0 * plane + i1 * n2 + 1; i2 < n2 - 1; i2++, v++) {
        const x = volume[v];
        if (!(x > level)) continue;
        let top = true;
        for (let k = 0; k < 26; k++) {
          const y = volume[v + offsets[k]];
          if (y > x || (y === x && offsets[k] < 0)) { top = false; break; }
        }
        if (top) { index.push(v); value.push(x); }
      }
    }
    onProgress((i0 + 1) / n0);
  }
  const order = index.map((_, k) => k).sort((p, q) => value[q] - value[p]);
  const halfWidth = halfWidths(volume, n, order.map((k) => index[k]), level0);
  const widest = Math.max(...halfWidth.filter(Number.isFinite));
  const R = radius ?? (Number.isFinite(widest) ? Math.min(12, Math.max(3, 1 + Math.ceil(2 * widest))) : 4);
  // A peak is kept when no stronger one is centred within 2R - 1 voxels of its maximum or of
  // its own centre (two maxima of one broad peak lead to the same centre).
  const peaks = [], taken = [];
  const crowded = (p) => taken.some((t) => Math.max(Math.abs(t[0] - p[0]), Math.abs(t[1] - p[1]), Math.abs(t[2] - p[2])) <= 2 * R - 1);
  let edge = 0, near = 0;
  for (const k of order) {
    if (peaks.length >= maxPeaks) break;
    const v = index[k], j = [v % n2, Math.floor(v / n2) % n1, Math.floor(v / plane)];
    if (crowded(j)) { near++; continue; }
    const fit = centroid(volume, n, j, R);
    if (!fit) { edge++; continue; }
    if (crowded(fit.c)) { near++; continue; }
    taken.push(fit.c);
    const u = fit.c.map((x, d) => axes[d].e0 + (x + 0.5) * axes[d].step);
    peaks.push({ h: mulv(W, u), u, max: value[k], intensity: fit.intensity, background: fit.background });
  }
  return { peaks, stats: { level: level0, spread, threshold: level, maxima: index.length, halfWidth, radius: R, edge, near } };
}

/** Each Q (Å⁻¹, no 2π) indexed with `ub`: { h, hkl (nearest integers), ok: within `tol` of them and not 000 }. */
export function indexQ(ub, qs, tol) {
  const inv = inv3(ub);
  return qs.map((q) => {
    const h = mulv(inv, q), hkl = h.map((x) => Math.round(x) || 0);
    const dev = Math.max(...h.map((x, i) => Math.abs(x - hkl[i])));
    return { h, hkl, dev, ok: dev <= tol && hkl.some(Boolean) };
  });
}

/** The volume of a cell (Å³). */
export function cellVolume([a, b, c, al, be, ga]) {
  const [ca, cb, cg] = [al, be, ga].map((x) => Math.cos(x * DEG));
  return a * b * c * Math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg);
}

/** Least squares over the indexed peaks: Q′ = R(p0..p2)·U0·B(cell(p3..)) with residuals in units of `scale`. */
function fitPeaks(qs, hkls, U0, sys, start, scale) {
  const model = (p) => mul(mul(smallRot([p[0], p[1], p[2]]), U0), bMatrix(...sys.cell(Array.from(p.slice(3)))));
  const resid = (p) => {
    const ub = model(p), r = new Float64Array(3 * hkls.length);
    hkls.forEach((h, k) => {
      const q = mulv(ub, h);
      for (let i = 0; i < 3; i++) r[3 * k + i] = (q[i] - qs[k][i]) / scale;
    });
    return r;
  };
  const steps = [1e-3, 1e-3, 1e-3, ...sys.names.map((name, i) => (name.length > 1 ? 1e-3 : start[i] * 1e-5))];
  const res = levenbergMarquardt(resid, [0, 0, 0, ...start], { steps, fScale: 1 });
  return { ub: model(res.p), p: res.p, covariance: res.covariance };
}

/**
 * Refine `ub0` (row-major, no 2π) against `peaks` (from findPeaks, at indices h of ub0) with
 * the cell constrained to `system`: index the peaks within `tol` of integers, fit, and
 * index again with the result, leaving out peaks whose residual |UB′·h − Q| exceeds 5 times
 * the median one (at least half of `scale`), until the set of peaks used stays (at most
 * `rounds` fits). Residuals are in units of `scale` (Å⁻¹, about a voxel), beyond which the
 * soft-L1 loss lowers a peak's weight. Returns { system, ub, cell, su (standard
 * uncertainties, 0 for fixed values), volume, volumeSU, rotation (°, the orientation's
 * change), rms (Å⁻¹), rmsH (r.l.u. per axis), indexed (one entry per peak: hkl, h under the
 * refined UB, dq in Å⁻¹, ok: within `tol`, used), within (peaks within `tol`), used,
 * cut (Å⁻¹), rounds, settled }.
 */
export function refineUB(peaks, ub0, system, { tol = 0.15, scale = 1e-3, rounds = 6 } = {}) {
  const sys = UB_SYSTEMS[system];
  if (!sys) throw new Error(`Unknown crystal system "${system}".`);
  if (!(det3(ub0) > 0)) throw new Error('The UB matrix is left-handed, so it has no orientation to refine.');
  const cell0 = cellFromUB(ub0);
  if (!fitsSystem(system, cell0, LOOSE)) {
    throw new Error(`The cell (a ${cell0[0].toFixed(4)}, b ${cell0[1].toFixed(4)}, c ${cell0[2].toFixed(4)} Å, ${cell0.slice(3).map((x) => x.toFixed(2)).join('/')}°) `
      + `does not have the metric of the ${sys.short} system in its standard setting; reindex it to that setting first.`);
  }
  const U0 = polarU(ub0, cell0, 1), start = startParams(system, cell0);
  const qs = peaks.map((p) => mulv(ub0, p.h));
  const miss = (u, k, hkl) => {
    const q = mulv(u, hkl);
    return Math.hypot(q[0] - qs[k][0], q[1] - qs[k][1], q[2] - qs[k][2]);
  };
  let ub = ub0, key = '', fit = null, set = [], hkls = [], cut = Infinity, round = 0, settled = false;
  for (; round < rounds; round++) {
    const idx = indexQ(ub, qs, tol);
    let next = idx.flatMap((x, k) => (x.ok ? [k] : []));
    if (next.length < 6) throw new Error(`Only ${next.length} of ${peaks.length} peaks index within ${tol} of integers; at least 6 are needed.`);
    // After a fit, peaks far off the others (a twin, a powder ring, a spurion) are left out.
    if (fit) {
      const dq = next.map((k) => miss(ub, k, idx[k].hkl));
      cut = Math.max(5 * median(dq), 0.5 * scale);
      next = next.filter((_, i) => dq[i] <= cut);
    }
    const nextKey = next.map((k) => `${k}:${idx[k].hkl}`).join(' ');
    if (nextKey === key) { settled = true; break; }
    const hk = next.map((k) => idx[k].hkl);
    const M = [0, 1, 2].flatMap((i) => [0, 1, 2].map((j) => hk.reduce((s, h) => s + h[i] * h[j], 0)));
    if (!(det3(M) > 1e-6 * ((M[0] + M[4] + M[8]) / 3) ** 3)) throw new Error('The indexed peaks lie in one plane of reciprocal space, which fixes no cell.');
    [set, hkls, key] = [next, hk, nextKey];
    fit = fitPeaks(set.map((k) => qs[k]), hkls, U0, sys, start, scale);
    ub = fit.ub;
  }
  const p = Array.from(fit.p), np = p.length, C = fit.covariance;
  const cellOf = (x) => sys.cell(x.slice(3));
  const cell = cellOf(p);
  // σ of a function of the parameters, from its numerical gradient and their covariance.
  const su = (f) => {
    const g = p.map((x, k) => {
      const h = k < 3 ? 1e-4 : Math.max(1e-6, Math.abs(x) * 1e-6), q = p.slice();
      q[k] += h;
      return (f(q) - f(p)) / h;
    });
    let s = 0;
    for (let a = 0; a < np; a++) for (let b = 0; b < np; b++) s += g[a] * C[a * np + b] * g[b];
    return Math.sqrt(Math.max(0, s));
  };
  const R = smallRot(p.slice(0, 3)), inSet = new Set(set);
  const indexed = indexQ(ub, qs, tol).map((x, k) => ({ hkl: x.hkl, h: x.h, dq: miss(ub, k, x.hkl), ok: x.ok, used: inSet.has(k) }));
  const used = indexed.filter((x) => x.used);
  return {
    system, ub, cell, su: cell.map((_, i) => su((q) => cellOf(q)[i])), volume: cellVolume(cell), volumeSU: su((q) => cellVolume(cellOf(q))),
    rotation: Math.acos(Math.max(-1, Math.min(1, (R[0] + R[4] + R[8] - 1) / 2))) / DEG,
    rms: Math.sqrt(used.reduce((s, x) => s + x.dq * x.dq, 0) / used.length),
    rmsH: [0, 1, 2].map((i) => Math.sqrt(used.reduce((s, x) => s + (x.h[i] - x.hkl[i]) ** 2, 0) / used.length)),
    indexed, within: indexed.filter((x) => x.ok).length, used: used.length, cut, rounds: round, settled,
  };
}

/** Q (Å⁻¹, no 2π) of one voxel step along each display axis of `dims` under `ub`: the scale of the residuals. */
export function voxelQ(dims, ub) {
  const W = basisMatrix(dims), axes = uniformAxes(dims);
  return axes.map((a, d) => {
    const q = mulv(ub, [W[d], W[3 + d], W[6 + d]]);
    return Math.hypot(...q) * a.step;
  });
}

/**
 * `x` with its standard uncertainty in parentheses, as crystallographers write it: the
 * s.u. to two digits when they are 19 or less, else to one (4.20003(2), 6.8001(15)); `x`
 * to `digits` decimals when `su` is 0 or unknown.
 */
export function formatSU(x, su, digits = 4) {
  if (!(su > 0) || !Number.isFinite(su)) return x.toFixed(digits);
  let d = 1 - Math.floor(Math.log10(su)), n = Math.round(su * 10 ** d);
  if (n > 19) { d -= 1; n = Math.round(su * 10 ** d); }
  return d >= 0 ? `${x.toFixed(d)}(${n})` : `${Math.round(x)}(${Math.round(su)})`;
}

/**
 * The peaks of a refinement as CSV, one row per peak found: its integer indices, its indices
 * in the volume (under the volume's UB) and under the refined UB, d (Å), the residual
 * |UB′·h − Q| (Å⁻¹, no 2π), its maximum and windowed intensity above background, and
 * whether it was used.
 */
export function peaksCSV(peaks, result) {
  const rows = [['H', 'K', 'L', 'h_volume', 'k_volume', 'l_volume', 'h_refined', 'k_refined', 'l_refined', 'd_A', 'residual_invA_no2pi', 'max', 'intensity', 'used'].join(',')];
  const num = (x, p) => (Number.isFinite(x) ? Number(x.toPrecision(p)).toString() : '');
  peaks.forEach((p, k) => {
    const x = result.indexed[k], q = mulv(result.ub, x.h);
    rows.push([...x.hkl, ...p.h.map((v) => num(v, 7)), ...x.h.map((v) => num(v, 7)), num(1 / Math.hypot(...q), 6), num(x.dq, 4),
      num(p.max, 6), num(p.intensity, 6), x.used ? 1 : 0].join(','));
  });
  return `${rows.join('\n')}\n`;
}

/** The largest change |T·h - h| of any index over the corners of the grid of `dims` (h′ = T·h). */
export function indexShift(T, dims) {
  const W = basisMatrix(dims), ends = dims.map((d) => [d.edges[0], d.edges.at(-1)]);
  let max = 0;
  for (let c = 0; c < 8; c++) {
    const h = mulv(W, [0, 1, 2].map((d) => ends[d][(c >> d) & 1]));
    const t = mulv(T, h);
    for (let i = 0; i < 3; i++) max = Math.max(max, Math.abs(t[i] - h[i]));
  }
  return max;
}
