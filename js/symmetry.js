// Point-group operations acting on reciprocal-space coordinates.
//
// An operation is a 3x3 integer matrix R (row-major) that maps a column
// (h, k, l) to R (h, k, l); it is written "h+k,-h,l" in the same way as the
// coordinate triplets of the International Tables. Real-space triplets in
// x, y, z are accepted too and converted with R = W^T, since the reflection
// h is equivalent to W^T h when x -> W x is a symmetry of the structure.
//
// For averaging, every operation of the closed group is turned into an
// integer affine map on bin indices, so symmetry-equivalent voxels are found
// exactly, without interpolation.

import { reciprocalMetric } from './nexus.js';

const MAX_ORDER = 48;
const IDENTITY = [1, 0, 0, 0, 1, 0, 0, 0, 1];

// Laue classes with generators in hkl form (hexagonal axes for trigonal).
export const PRESETS = [
  ['1', 'None'],
  ['-1', '-h,-k,-l'],
  ['2/m (b unique)', '-h,k,-l; -h,-k,-l'],
  ['2/m (c unique)', '-h,-k,l; -h,-k,-l'],
  ['mmm', '-h,-k,l; -h,k,-l; -h,-k,-l'],
  ['4/m', 'k,-h,l; -h,-k,-l'],
  ['4/mmm', 'k,-h,l; h,-k,-l; -h,-k,-l'],
  ['-3', 'k,-h-k,l; -h,-k,-l'],
  ['-3m1', 'k,-h-k,l; h,-h-k,-l; -h,-k,-l'],
  ['-31m', 'k,-h-k,l; -k,-h,-l; -h,-k,-l'],
  ['6/m', 'h+k,-h,l; -h,-k,-l'],
  ['6/mmm', 'h+k,-h,l; -k,-h,-l; -h,-k,-l'],
  ['m-3', 'k,l,h; -h,-k,l; -h,k,-l; -h,-k,-l'],
  ['m-3m', 'k,l,h; k,-h,l; -h,-k,-l'],
];

/** Parse "h+k,-h,l" (or "x-y,x,z") into a row-major integer matrix acting on hkl. */
export function parseOp(text) {
  const src = text.replace(/\s+/g, '').toLowerCase();
  const parts = src.split(',');
  if (parts.length !== 3 || parts.some((p) => !p)) throw new Error(`"${text}" needs three comma-separated components.`);
  const real = /[xyz]/.test(src);
  if (real && /[hkl]/.test(src)) throw new Error(`"${text}" mixes h,k,l with x,y,z.`);
  const vars = real ? 'xyz' : 'hkl';
  const m = [];
  for (const part of parts) {
    const row = [0, 0, 0];
    const term = /([+-]?)(\d*)\*?([a-z])/gy;
    let used = 0, match;
    while ((match = term.exec(part))) {
      const [, sign, num, v] = match;
      const i = vars.indexOf(v);
      if (i < 0) throw new Error(`Unknown symbol "${v}" in "${text}".`);
      row[i] += (sign === '-' ? -1 : 1) * (num === '' ? 1 : Number(num));
      used = term.lastIndex;
    }
    if (used !== part.length) {
      throw new Error(`Cannot read "${part}" in "${text}". Use integer combinations such as h-k or -2l; translations are not allowed.`);
    }
    m.push(...row);
  }
  const R = real ? transpose(m) : m;
  if (Math.abs(det(R)) !== 1) throw new Error(`"${text}" has determinant ${det(R)}; a point operation must have determinant ±1.`);
  return R;
}

/** Parse a list of operations separated by ";" or new lines. Empty text means identity only. */
export function parseOps(text) {
  return text.split(/[;\n]/).map((s) => s.trim()).filter((s) => s && s.toLowerCase() !== 'none').map(parseOp);
}

/** All products of the generators: the smallest group containing them. */
export function closeGroup(generators) {
  const ops = [IDENTITY], keys = new Set([IDENTITY.join()]);
  const add = (m) => {
    const key = m.join();
    if (keys.has(key)) return false;
    if (ops.length >= MAX_ORDER) {
      throw new Error(`These operations generate more than ${MAX_ORDER} operations, so they are not a crystallographic point group.`);
    }
    keys.add(key);
    ops.push(m);
    return true;
  };
  generators.forEach(add);
  for (let grew = true; grew;) {
    grew = false;
    for (const a of [...ops]) for (const b of [...ops]) grew = add(multiply(a, b)) || grew;
  }
  return ops;
}

/** "h+k,-h,l" */
export function formatOp(m) {
  const rows = [0, 1, 2].map((i) => {
    let s = '';
    for (let j = 0; j < 3; j++) {
      const c = m[3 * i + j];
      if (!c) continue;
      const mag = Math.abs(c) === 1 ? '' : String(Math.abs(c));
      s += (c < 0 ? '-' : s ? '+' : '') + mag + 'hkl'[j];
    }
    return s || '0';
  });
  return rows.join(',');
}

/**
 * Largest change of the reciprocal metric under any operation, relative to
 * its largest element. Zero for true symmetries of the cell.
 */
export function metricChange(ops, cell) {
  const G = reciprocalMetric(cell);
  const scale = Math.max(...G.flat().map(Math.abs));
  let worst = 0;
  for (const R of ops) {
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let v = 0;
        for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) v += R[3 * a + i] * G[a][b] * R[3 * b + j];
        worst = Math.max(worst, Math.abs(v - G[i][j]) / scale);
      }
    }
  }
  return worst;
}

/**
 * Express each operation as an integer affine map on display-dimension bin
 * indices: i' = M i + t. dims[d].basis.vec gives the HKL basis vector of
 * display axis d; without bases the operations act on the axes directly.
 * Throws when an operation does not send bin centers onto bin centers.
 */
export function indexMaps(ops, dims) {
  const B = dims.every((d) => d.basis) ? transpose(dims.flatMap((d) => d.basis.vec)) : IDENTITY;
  const Binv = inverse(B);
  const n = dims.map((d) => d.edges.length - 1);
  const w = dims.map((d, i) => (d.edges[n[i]] - d.edges[0]) / n[i]);
  const c0 = dims.map((d, i) => d.edges[0] + w[i] / 2);
  dims.forEach((d, i) => {
    for (let k = 0; k < n[i]; k++) {
      if (Math.abs(d.edges[k + 1] - d.edges[k] - w[i]) > 1e-4 * w[i]) {
        throw new Error(`Axis ${d.label} has non-uniform bins; symmetry averaging needs uniform bins.`);
      }
    }
  });
  const near = (x) => Math.abs(x - Math.round(x)) < 1e-3;
  return ops.map((R) => {
    const A = multiply(Binv, multiply(R, B));
    const M = new Int32Array(9), t = new Int32Array(3);
    for (let d = 0; d < 3; d++) {
      let offset = -c0[d];
      for (let e = 0; e < 3; e++) {
        const m = A[3 * d + e] * w[e] / w[d];
        if (!near(m)) {
          throw new Error(`Operation ${formatOp(R)} does not map the bin grid onto itself: axes ${dims[d].label} and ${dims[e].label} need compatible bin widths.`);
        }
        M[3 * d + e] = Math.round(m);
        offset += A[3 * d + e] * c0[e];
      }
      if (!near(offset / w[d])) {
        throw new Error(`Operation ${formatOp(R)} does not map bin centers onto bin centers along ${dims[d].label}: the grid needs a bin centered at 0.`);
      }
      t[d] = Math.round(offset / w[d]);
    }
    return { M, t };
  });
}

// ---- 3x3 helpers (row-major arrays of 9) -------------------------------------

function multiply(a, b) {
  const out = new Array(9).fill(0);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) out[3 * i + j] += a[3 * i + k] * b[3 * k + j];
  return out;
}

function transpose(m) {
  return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}

function det(m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

function inverse(m) {
  const d = det(m);
  if (Math.abs(d) < 1e-12) throw new Error('Axis basis vectors are linearly dependent.');
  const c = (i, j) => {
    const r = [0, 1, 2].filter((x) => x !== i), s = [0, 1, 2].filter((x) => x !== j);
    return ((i + j) % 2 ? -1 : 1) * (m[3 * r[0] + s[0]] * m[3 * r[1] + s[1]] - m[3 * r[0] + s[1]] * m[3 * r[1] + s[0]]);
  };
  const out = new Array(9);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) out[3 * i + j] = c(j, i) / d;
  return out;
}
