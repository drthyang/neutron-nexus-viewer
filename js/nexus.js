// Locate and read a 3-D histogram from a NeXus/HDF5 file opened with h5wasm.
//
// Primary target: Mantid MDHistoWorkspace files written by SaveMD
// (/MDHistoWorkspace/data/{signal,mask,D0,D1,D2}). Any NXdata group with a
// 3-D signal and bin-edge or bin-center axes also works. Size-1 dimensions
// are dropped, so a 4-D workspace with one integrated axis is accepted.
//
// Display dimension d lives on storage axis keep[2 - d]; for Mantid files
// (signal axes "D2:D1:D0") this makes dims[0] = D0, dims[1] = D1, dims[2] = D2.

const MAX_VOLUME_BYTES = 3.5e9;

const str = (v) => {
  if (v == null) return '';
  if (Array.isArray(v) || ArrayBuffer.isView(v)) return v.length ? str(v[0]) : '';
  return String(v).trim();
};
const list = (v) => {
  if (v == null) return null;
  if (Array.isArray(v)) return v.map(str);
  const s = str(v);
  return s ? s.split(/[:,]/).map((x) => x.trim()) : null;
};
const attr = (node, name) => node.attrs[name]?.value;
const isGroup = (node) => node?.type === 'Group';
const isDataset = (node) => node?.type === 'Dataset';
const join = (base, name) => (base === '/' ? '' : base) + '/' + name;
const scalar = (node) => Number(isDataset(node) ? toNumbers(node.value)[0] : NaN);

function toNumbers(v) {
  if (v instanceof BigInt64Array || v instanceof BigUint64Array) return Float64Array.from(v, Number);
  if (ArrayBuffer.isView(v) || Array.isArray(v)) return v;
  return [Number(v)];
}

function* groups(node, path, depth) {
  yield [node, path];
  if (depth === 0) return;
  for (const key of node.keys()) {
    let child;
    try { child = node.get(key); } catch { continue; }
    if (isGroup(child)) yield* groups(child, join(path, key), depth - 1);
  }
}

/** Describe the largest 3-D histogram in the file, or throw a readable error. */
export function describeFile(file) {
  const found = [];
  const problems = [];
  for (const [group, path] of groups(file, '/', 4)) {
    if (str(attr(group, 'NX_class')) !== 'NXdata') continue;
    try {
      const info = describeNXdata(group, path);
      if (info) found.push(info);
    } catch (err) {
      problems.push(`${path}: ${err.message}`);
    }
  }
  if (!found.length) {
    throw new Error(problems.length ? problems.join(' ')
      : 'No NXdata group with a 3-D signal was found. Expected a Mantid MDHistoWorkspace saved with SaveMD.');
  }
  const info = found.sort((a, b) => b.voxels - a.voxels)[0];
  const entry = info.group.split('/').slice(0, -1).join('/') || '/';
  info.lattice = findLattice(file, entry);
  return info;
}

function describeNXdata(group, path) {
  const keys = group.keys();
  let name = str(attr(group, 'signal'));
  if (!keys.includes(name)) {
    name = keys.find((k) => {
      const node = group.get(k);
      return isDataset(node) && str(attr(node, 'signal')) === '1';
    });
  }
  if (!name) return null;
  const signal = group.get(name);
  const shape = signal.shape ?? [];
  const keep = shape.map((n, i) => (n > 1 ? i : -1)).filter((i) => i >= 0);
  if (keep.length !== 3) {
    throw new Error(`signal "${name}" has ${keep.length} non-trivial dimensions; this viewer needs 3.`);
  }
  // HDF5 type class 0 is integer, 1 is float.
  if (![0, 1].includes(signal.metadata?.type)) {
    throw new Error(`signal "${name}" is not numeric (${JSON.stringify(signal.dtype)}).`);
  }
  const voxels = keep.reduce((n, a) => n * shape[a], 1);
  if (voxels * 4 > MAX_VOLUME_BYTES) {
    throw new Error(`signal is ${(voxels * 4 / 1e9).toFixed(1)} GB as float32, too large to hold in browser memory. Rebin it coarser first.`);
  }

  const axisNames = list(attr(group, 'axes')) ?? list(attr(signal, 'axes'));
  const maskNode = keys.includes('mask') ? group.get('mask') : null;
  const mask = isDataset(maskNode) && String(maskNode.shape) === String(shape) ? join(path, 'mask') : null;

  const dims = [2, 1, 0].map((k) => {
    const axis = keep[k], n = shape[axis];
    const axisName = axisNames?.length === shape.length ? axisNames[axis] : null;
    const node = axisName && keys.includes(axisName) ? group.get(axisName) : null;
    let edges = Array.from({ length: n + 1 }, (_, i) => i);
    let longName = axisName && axisName !== '.' ? axisName : `axis ${axis}`, units = '', frame = '';
    if (isDataset(node)) {
      const values = Array.from(toNumbers(node.value), Number);
      if (values.length === n + 1) edges = values;
      else if (values.length === n) edges = centersToEdges(values);
      else throw new Error(`axis "${axisName}" has ${values.length} values for ${n} bins.`);
      longName = str(attr(node, 'long_name')) || axisName;
      units = str(attr(node, 'units'));
      frame = str(attr(node, 'frame'));
    }
    if (edges[n] < edges[0]) throw new Error(`axis "${longName}" is descending; only ascending axes are supported.`);
    const basis = parseBasis(longName);
    return { name: axisName ?? `axis${axis}`, longName, units, frame, edges, basis, label: shortLabel(longName, basis) };
  });

  const chunks = signal.metadata?.chunks ?? null;
  return {
    group: path,
    signal: join(path, name),
    mask,
    dtype: signal.dtype,
    storageShape: shape,
    keep,
    shape: keep.map((a) => shape[a]),
    chunks,
    filters: (signal.filters ?? []).map((f) => f.name),
    voxels,
    dims,
  };
}

function centersToEdges(c) {
  const n = c.length;
  if (n === 1) return [c[0] - 0.5, c[0] + 0.5];
  const e = [c[0] - (c[1] - c[0]) / 2];
  for (let i = 0; i + 1 < n; i++) e.push((c[i] + c[i + 1]) / 2);
  e.push(c[n - 1] + (c[n - 1] - c[n - 2]) / 2);
  return e;
}

/** "[H,0,0]" -> {vec: [1,0,0], letter: "H"}; "[-K,K,0]" -> {vec: [-1,1,0], letter: "K"}. */
export function parseBasis(longName) {
  const m = /^\[(.*)\]$/.exec(longName.trim());
  if (!m) return null;
  const parts = m[1].split(',').map((s) => s.replace(/\s+/g, ''));
  if (parts.length !== 3) return null;
  const vec = [], letters = new Set();
  for (const part of parts) {
    const q = /^([+-]?(?:\d+(?:\.\d*)?|\.\d+)?)\*?([A-Za-z]*)$/.exec(part);
    if (!q) return null;
    const [, num, letter] = q;
    if (!letter) {
      if (num === '' || Number(num) !== 0) return null;
      vec.push(0);
    } else {
      letters.add(letter.toUpperCase());
      vec.push(num === '' || num === '+' ? 1 : num === '-' ? -1 : Number(num));
    }
  }
  if (letters.size !== 1 || vec.every((v) => v === 0)) return null;
  return { vec, letter: [...letters][0] };
}

function shortLabel(longName, basis) {
  if (!basis) return longName;
  const i = 'HKL'.indexOf(basis.letter);
  const unit = i >= 0 && basis.vec.every((v, j) => v === (j === i ? 1 : 0));
  return unit ? basis.letter : longName;
}

/**
 * Unit cell of the sample. Prefers the UB matrix (Mantid oriented_lattice
 * orientation_matrix), whose metric (UB)^T UB is the reciprocal metric without
 * 2*pi; falls back to stored unit_cell values.
 */
function findLattice(file, entryPath) {
  const entry = file.get(entryPath);
  if (!isGroup(entry)) return null;
  let stored = null;
  for (const [group] of groups(entry, entryPath, 4)) {
    const keys = group.keys();
    if (keys.includes('orientation_matrix')) {
      const ub = Array.from(toNumbers(group.get('orientation_matrix').value ?? []), Number);
      const cell = ub.length === 9 && ub.every(Number.isFinite) ? cellFromUB(ub) : null;
      if (cell) return { ...cell, source: 'UB', ub };
    }
    let values = null;
    if (keys.includes('unit_cell_a')) {
      values = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].map((k) => scalar(group.get(`unit_cell_${k}`)));
    } else if (keys.includes('unit_cell')) {
      const v = Array.from(toNumbers(group.get('unit_cell').value ?? []), Number);
      if (v.length === 6) values = v;
    }
    if (!stored && values && values.every(Number.isFinite) && values.slice(0, 3).every((x) => x > 0)
        && values.slice(3).every((x) => x > 0 && x < 180)) {
      const [a, b, c, alpha, beta, gamma] = values;
      stored = { a, b, c, alpha, beta, gamma, source: 'unit cell' };
    }
  }
  return stored;
}

/** Direct cell (angstrom, degrees) from a row-major UB matrix; null if singular. */
export function cellFromUB(ub) {
  const gstar = [0, 1, 2].map((i) => [0, 1, 2].map((j) => ub[i] * ub[j] + ub[3 + i] * ub[3 + j] + ub[6 + i] * ub[6 + j]));
  const g = invert3(gstar);
  if (!g) return null;
  const [a, b, c] = [0, 1, 2].map((i) => Math.sqrt(g[i][i]));
  const angle = (i, j) => Math.acos(g[i][j] / Math.sqrt(g[i][i] * g[j][j])) * 180 / Math.PI;
  const cell = { a, b, c, alpha: angle(1, 2), beta: angle(0, 2), gamma: angle(0, 1) };
  return Object.values(cell).every(Number.isFinite) ? cell : null;
}

function invert3(m) {
  const det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  if (!(Math.abs(det) > 1e-300)) return null;
  const inv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const [i1, i2] = [(j + 1) % 3, (j + 2) % 3], [j1, j2] = [(i + 1) % 3, (i + 2) % 3];
      inv[i][j] = (m[i1][j1] * m[i2][j2] - m[i1][j2] * m[i2][j1]) / det;
    }
  }
  return inv;
}

/**
 * Read signal (and mask) into a float32 volume in storage order, with masked
 * and non-finite voxels set to NaN. Reads whole chunks along the leading axis.
 */
export function loadVolume(file, info, onProgress = () => {}) {
  const signal = file.get(info.signal);
  const mask = info.mask ? file.get(info.mask) : null;
  const lead = info.keep[0];
  const [n0, n1, n2] = info.shape;
  const plane = n1 * n2;
  const chunk = info.chunks?.[lead] ?? 0;
  const step = Math.min(n0, Math.max(1, chunk || Math.floor(32e6 / (plane * 8))));
  let volume;
  try {
    volume = new Float32Array(n0 * plane);
  } catch {
    throw new Error(`Could not allocate ${(n0 * plane * 4 / 1e9).toFixed(2)} GB for the volume.`);
  }
  let valid = 0, min = Infinity, max = -Infinity;
  for (let i = 0; i < n0; i += step) {
    const j = Math.min(n0, i + step);
    const ranges = info.storageShape.map((n, a) => (a === lead ? [i, j] : n === 1 ? [0, 1] : []));
    const s = toNumbers(signal.slice(ranges));
    const m = mask ? toNumbers(mask.slice(ranges)) : null;
    const offset = i * plane;
    for (let k = 0; k < s.length; k++) {
      const v = s[k];
      if ((m && m[k] != 0) || !Number.isFinite(v)) {
        volume[offset + k] = NaN;
      } else {
        volume[offset + k] = v;
        valid++;
        if (v < min) min = v;
        if (v > max) max = v;
      }
    }
    onProgress(j / n0);
  }
  return { volume, stats: { valid, fraction: valid / volume.length, min, max } };
}

// ---- Geometry -------------------------------------------------------------

const NOMINAL_ANGLES = [60, 90, 120];

/** Snap direct-cell angles within 1 degree of 60/90/120 to those values. */
export function nominalCell(cell) {
  const snap = (x) => NOMINAL_ANGLES.find((t) => Math.abs(x - t) <= 1) ?? x;
  return { ...cell, alpha: snap(cell.alpha), beta: snap(cell.beta), gamma: snap(cell.gamma) };
}

/** Reciprocal metric tensor (no 2*pi) of a direct cell given in angstrom and degrees. */
export function reciprocalMetric({ a, b, c, alpha, beta, gamma }) {
  const r = Math.PI / 180;
  const ca = Math.cos(alpha * r), cb = Math.cos(beta * r), cg = Math.cos(gamma * r);
  return invert3([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]]);
}

/**
 * Display->Cartesian matrix (row-major 3x3) for the 3-D view. With a cell and
 * HKL bases it applies the reciprocal lattice vectors (no 2*pi), with c* along
 * z and a* in the x-z plane; otherwise axes sharing a unit keep unit scale and
 * mixed units are stretched to a common extent.
 */
export function cartesianBasis(dims, cell) {
  if (cell && dims.every((d) => d.basis && isHKL(d))) {
    const G = reciprocalMetric(cell);
    const [as, bs, cs] = [0, 1, 2].map((i) => Math.sqrt(G[i][i]));
    const cosAlpha = G[1][2] / (bs * cs), cosBeta = G[0][2] / (as * cs);
    const a = [as * Math.sqrt(1 - cosBeta * cosBeta), 0, as * cosBeta];
    const bz = bs * cosAlpha, bx = (G[0][1] - a[2] * bz) / a[0];
    const L = [a, [bx, Math.sqrt(Math.max(0, bs * bs - bx * bx - bz * bz)), bz], [0, 0, cs]];
    const T = [];
    for (let i = 0; i < 3; i++) for (let d = 0; d < 3; d++) T.push(dims[d].basis.vec.reduce((s, h, k) => s + h * L[k][i], 0));
    return { T, lattice: true };
  }
  const extent = dims.map((d) => d.edges[d.edges.length - 1] - d.edges[0]);
  const shared = dims[0].units && dims.every((d) => d.units === dims[0].units);
  const s = shared ? [1, 1, 1] : extent.map((e) => Math.max(...extent) / e);
  return { T: [s[0], 0, 0, 0, s[1], 0, 0, 0, s[2]], lattice: false };
}

const isHKL = (dim) => dim.frame === 'HKL' || (!dim.frame && /r\.?l\.?u/i.test(dim.units));

/**
 * Display geometry of the plane spanned by dims x and y: axis lengths per unit
 * coordinate, the cosine of the angle between them, and whether both axes
 * share a unit (equal aspect). Falls back to rectangular axes without a cell.
 */
export function planeGeometry(dims, cell, x, y) {
  const X = dims[x], Y = dims[y];
  if (cell && X.basis && Y.basis && isHKL(X) && isHKL(Y)) {
    const G = reciprocalMetric(cell);
    const dot = (u, v) => u.reduce((s, ui, i) => s + ui * v.reduce((t, vj, j) => t + G[i][j] * vj, 0), 0);
    const lx = Math.sqrt(dot(X.basis.vec, X.basis.vec)), ly = Math.sqrt(dot(Y.basis.vec, Y.basis.vec));
    return { lx, ly, cos: dot(X.basis.vec, Y.basis.vec) / (lx * ly), equal: true, lattice: true };
  }
  return { lx: 1, ly: 1, cos: 0, equal: !!X.units && X.units === Y.units, lattice: false };
}
