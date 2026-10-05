// Locate and read a 3-D histogram from a NeXus/HDF5 file opened with h5wasm.
//
// Primary target: Mantid MDHistoWorkspace files written by SaveMD
// (/MDHistoWorkspace/data/{signal,mask,D0,D1,D2}). Any NXdata group with a
// 3-D signal and bin-edge or bin-center axes also works. Size-1 dimensions
// are dropped, so a 4-D workspace with one integrated axis is accepted.
// Files without NXdata written by NEBULA3D (https://github.com/drthyang/nebula3d)
// open too: its reciprocal-space volumes (/entry/{data, sigma, mask, h_axis,
// k_axis, l_axis, ub_matrix}) and its real-space 3D-ΔPDFs (/{data, x_axis,
// y_axis, z_axis} with the cell in lat_* attributes).
//
// The volume is held in the order of info.shape, the storage axes listed in
// info.keep, and display dimension d lives on storage axis keep[2 - d]. For
// Mantid files (signal axes "D2:D1:D0") keep is in storage order, which makes
// dims[0] = D0, dims[1] = D1, dims[2] = D2. When keep is not in storage order
// (a dataset aligned to another's axes, see alignDims()), the volume is a
// transpose of the stored array, reordered as it is read.

const MAX_VOLUME_BYTES = 3.5e9;
const NEBULA3D_AXES = [
  { names: ['h_axis', 'k_axis', 'l_axis'], longNames: ['[H,0,0]', '[0,K,0]', '[0,0,L]'], units: 'r.l.u.', frame: 'HKL' },
  { names: ['x_axis', 'y_axis', 'z_axis'], longNames: ['x', 'y', 'z'], units: 'Å', frame: 'direct' },
];
// Real-space axes along a, b and c in Å, as NEBULA3D writes its 3D-ΔPDFs in the Mantid layout.
const DIRECT_NAME = /^[xyz]$/i, ANGSTROM = /^(Å|Angstroms?)$/i;

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

/**
 * Describe the largest 3-D histogram in the file, or throw a readable error.
 * With `order` (the labels of another dataset's display axes), axes that are
 * the same set in another order are reordered to match it.
 */
export function describeFile(file, { order = null } = {}) {
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
  if (!found.length && !problems.length) {
    const info = describeNebula(file);
    if (info) return order ? alignDims(info, order) : info;
  }
  if (!found.length) {
    throw new Error(problems.length ? problems.join(' ')
      : 'No 3-D histogram was found: expected a Mantid MDHistoWorkspace saved with SaveMD, an NXdata group with a 3-D signal, or a NEBULA3D volume or 3D-ΔPDF.');
  }
  const best = found.sort((a, b) => b.voxels - a.voxels)[0];
  const entry = best.group.split('/').slice(0, -1).join('/') || '/';
  const info = withCell(best, findLattice(file, entry));
  return order ? alignDims(info, order) : info;
}

/**
 * Attach the unit cell (or null): it gives axes along the direct lattice their
 * length in Å, and marks real-space data `signed` (a ΔPDF changes sign, so it
 * is shown on a range symmetric about 0).
 */
function withCell(info, lattice) {
  const real = info.dims.every(isDirect);
  const dims = real ? info.dims.map((d) => ({ ...d, length: lattice ? [lattice.a, lattice.b, lattice.c][d.axis] : null })) : info.dims;
  return { ...info, dims, lattice, signed: real };
}

/**
 * Display `info`'s axes in the order of the labels `order` when they are the
 * same three axes, by changing which storage axis each display dimension reads.
 */
export function alignDims(info, order) {
  const at = order.map((label) => info.dims.findIndex((d) => d.label === label));
  if (at.some((i) => i < 0) || new Set(at).size !== 3 || at.every((i, d) => i === d)) return info;
  const keep = [0, 1, 2].map((j) => info.keep[2 - at[2 - j]]);
  return { ...info, keep, shape: keep.map((a) => info.storageShape[a]), dims: at.map((i) => info.dims[i]) };
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
  const sameShape = (key) => {
    const node = keys.includes(key) ? group.get(key) : null;
    return isDataset(node) && String(node.shape) === String(shape);
  };
  // Mantid's mask: nonzero marks a masked voxel.
  const mask = sameShape('mask') ? { path: join(path, 'mask'), valid: false } : null;
  // Uncertainties: Mantid writes variances (errors_squared), NeXus standard deviations (errors).
  const squared = ['errors_squared', `${name}_errors_squared`].find(sameShape);
  const plain = ['errors', `${name}_errors`].find(sameShape);
  const errors = squared ? { path: join(path, squared), squared: true } : plain ? { path: join(path, plain), squared: false } : null;

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
  // Axes x, y and z in Å (a 3D-ΔPDF) lie along the direct lattice vectors a, b and c.
  const xyz = dims.map((d) => (DIRECT_NAME.test(d.longName) && ANGSTROM.test(d.units) ? 'xyz'.indexOf(d.longName.toLowerCase()) : -1));
  if (xyz.every((a) => a >= 0) && new Set(xyz).size === 3) {
    dims.forEach((d, i) => Object.assign(d, { frame: 'direct', units: 'Å', axis: xyz[i], label: 'xyz'[xyz[i]] }));
  }

  const chunks = signal.metadata?.chunks ?? null;
  return {
    group: path,
    signal: join(path, name),
    mask,
    errors,
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

/**
 * A NEBULA3D volume or 3D-ΔPDF: a group with a 3-D `data` and one bin-center
 * axis per storage axis, h_axis, k_axis, l_axis (r.l.u.) or x_axis, y_axis,
 * z_axis (Å along the direct axes a, b, c), in C order with H or x slowest.
 * Volumes add `mask` (nonzero = valid), `sigma` (standard deviations) and
 * `ub_matrix` (with 2π); ΔPDFs store the direct cell as lat_* attributes.
 * Returns null when the file has no such group.
 */
function describeNebula(file) {
  for (const [group, path] of groups(file, '/', 2)) {
    const keys = group.keys();
    const kind = NEBULA3D_AXES.find((k) => k.names.every((name) => keys.includes(name)));
    if (!kind || !keys.includes('data')) continue;
    const signal = group.get('data');
    const shape = signal?.shape ?? [];
    if (!isDataset(signal) || shape.length !== 3) continue;
    if (![0, 1].includes(signal.metadata?.type)) throw new Error(`${join(path, 'data')} is not numeric (${JSON.stringify(signal.dtype)}).`);
    const voxels = shape[0] * shape[1] * shape[2];
    if (voxels * 4 > MAX_VOLUME_BYTES) {
      throw new Error(`data is ${(voxels * 4 / 1e9).toFixed(1)} GB as float32, too large to hold in browser memory.`);
    }
    const sameShape = (key) => keys.includes(key) && isDataset(group.get(key)) && String(group.get(key).shape) === String(shape);
    const real = kind.frame === 'direct';
    // Storage axis a holds H (or x) for a = 0. As for Mantid's D0, D1, D2, the
    // fastest storage axis is displayed first: L, K, H (or z, y, x).
    const dims = [2, 1, 0].map((a) => {
      const name = kind.names[a], values = Array.from(toNumbers(group.get(name).value), Number);
      if (values.length !== shape[a]) throw new Error(`${join(path, name)} has ${values.length} values for ${shape[a]} bins.`);
      const edges = centersToEdges(values);
      if (edges[shape[a]] < edges[0]) throw new Error(`${join(path, name)} is descending; only ascending axes are supported.`);
      const longName = kind.longNames[a], basis = real ? null : parseBasis(longName);
      const dim = { name, longName, units: kind.units, frame: kind.frame, edges, basis, label: real ? longName : shortLabel(longName, basis) };
      return real ? { ...dim, axis: a } : dim;
    });
    return withCell({
      group: path,
      signal: join(path, 'data'),
      mask: sameShape('mask') ? { path: join(path, 'mask'), valid: true } : null,
      errors: sameShape('sigma') ? { path: join(path, 'sigma'), squared: false } : null,
      dtype: signal.dtype,
      storageShape: shape,
      keep: [0, 1, 2],
      shape,
      chunks: signal.metadata?.chunks ?? null,
      filters: (signal.filters ?? []).map((f) => f.name),
      voxels,
      dims,
    }, real ? cellFromAttrs(group) : nebulaLattice(group));
  }
  return null;
}

/** The cell of a NEBULA3D volume from its ub_matrix (with 2π); null without one, or for the identity it stores when the UB is unknown. */
function nebulaLattice(group) {
  if (!group.keys().includes('ub_matrix')) return null;
  const stored = Array.from(toNumbers(group.get('ub_matrix').value ?? []), Number);
  if (stored.length !== 9 || stored.every((x, i) => x === (i % 4 ? 0 : 1))) return null;
  const ub = stored.map((x) => x / (2 * Math.PI));
  const cell = ub.every(Number.isFinite) ? cellFromUB(ub) : null;
  return cell ? { ...cell, source: 'UB', ub } : null;
}

/** The direct cell in the lat_a, lat_b, lat_c (Å) and lat_alpha, lat_beta, lat_gamma (degrees) attributes; null unless all six are valid. */
function cellFromAttrs(group) {
  const values = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].map((k) => Number(toNumbers(attr(group, `lat_${k}`) ?? NaN)[0]));
  if (!values.every(Number.isFinite) || !values.slice(0, 3).every((x) => x > 0) || !values.slice(3).every((x) => x > 0 && x < 180)) return null;
  const [a, b, c, alpha, beta, gamma] = values;
  return { a, b, c, alpha, beta, gamma, source: 'lat_* attributes' };
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
 * Read signal (and mask) into a float32 volume in the order of info.shape,
 * with masked and non-finite voxels set to NaN. Reads whole chunks along the
 * leading storage axis.
 */
export function loadVolume(file, info, onProgress = () => {}) {
  const signal = file.get(info.signal);
  const mask = info.mask ? file.get(info.mask.path) : null;
  // Mantid marks masked voxels with nonzero values, NEBULA3D valid ones.
  const maskedIf = !info.mask?.valid;
  const volume = allocate(info, 'the volume');
  let valid = 0, min = Infinity, max = -Infinity;
  readChunks(info, volume, (ranges, out) => {
    const s = toNumbers(signal.slice(ranges));
    const m = mask ? toNumbers(mask.slice(ranges)) : null;
    for (let k = 0; k < s.length; k++) {
      const v = s[k];
      if ((m && (m[k] != 0) === maskedIf) || !Number.isFinite(v)) {
        out[k] = NaN;
      } else {
        out[k] = v;
        valid++;
        if (v < min) min = v;
        if (v > max) max = v;
      }
    }
  }, onProgress);
  return { volume, stats: { valid, fraction: valid / volume.length, min, max } };
}

/**
 * Per-voxel variances σ² in the order of info.shape, from the errors dataset
 * (squared when it holds standard deviations), with NaN where unknown; null
 * when the file has no uncertainties. Errors that are zero everywhere (as
 * Mantid's layout requires them, written for a 3D-ΔPDF) count as none.
 */
export function loadVariance(file, info, onProgress = () => {}) {
  if (!info.errors) return null;
  const errors = file.get(info.errors.path), squared = info.errors.squared;
  const variance = allocate(info, 'the uncertainties');
  let positive = false;
  readChunks(info, variance, (ranges, out) => {
    const e = toNumbers(errors.slice(ranges));
    for (let k = 0; k < e.length; k++) {
      const v = squared ? e[k] : e[k] * e[k];
      out[k] = v >= 0 && v < Infinity ? v : NaN;
      if (v > 0) positive = true;
    }
  }, onProgress);
  return positive ? variance : null;
}

function allocate(info, what) {
  const [n0, n1, n2] = info.shape;
  try {
    return new Float32Array(n0 * n1 * n2);
  } catch {
    throw new Error(`Could not allocate ${(n0 * n1 * n2 * 4 / 1e9).toFixed(2)} GB for ${what}.`);
  }
}

/**
 * Fill `target` (in the order of info.shape) from datasets shaped like the
 * signal, in whole chunks along the leading storage axis: fill(ranges, out)
 * writes the slab `ranges` into `out` in storage order, and readChunks puts it
 * in place, reordering its axes when info.keep is not in storage order.
 */
function readChunks(info, target, fill, onProgress) {
  const { keep, shape, storageShape } = info;
  const axes = [...keep].sort((a, b) => a - b);
  const lead = axes[0], n = storageShape[lead];
  const plane = target.length / n;
  const chunk = info.chunks?.[lead] ?? 0;
  const step = Math.min(n, Math.max(1, chunk || Math.floor(32e6 / (plane * 8))));
  // The stride in `target` of each storage axis, slowest first.
  const [s0, s1, s2] = axes.map((a) => shape.slice(keep.indexOf(a) + 1).reduce((p, m) => p * m, 1));
  const [, m1, m2] = axes.map((a) => storageShape[a]);
  const inOrder = s0 > s1 && s1 > s2;
  const slab = inOrder ? null : new Float32Array(step * plane);
  for (let i = 0; i < n; i += step) {
    const j = Math.min(n, i + step);
    const ranges = storageShape.map((m, a) => (a === lead ? [i, j] : m === 1 ? [0, 1] : []));
    if (inOrder) {
      fill(ranges, target.subarray(i * plane, j * plane));
    } else {
      fill(ranges, slab);
      for (let a = i, k = 0; a < j; a++) {
        for (let b = 0; b < m1; b++) {
          for (let c = 0, o = a * s0 + b * s1; c < m2; c++, o += s2) target[o] = slab[k++];
        }
      }
    }
    onProgress(j / n);
  }
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
  if (cell && dims.every(isDirect)) {
    // Coordinates in Å along a, b and c: their unit vectors, a along x and b in the x-y plane.
    const r = Math.PI / 180, [ca, cb, cg] = [cell.alpha, cell.beta, cell.gamma].map((x) => Math.cos(x * r));
    const sg = Math.sin(cell.gamma * r), cy = (ca - cb * cg) / sg;
    const L = [[1, 0, 0], [cg, sg, 0], [cb, cy, Math.sqrt(Math.max(0, 1 - cb * cb - cy * cy))]];
    const T = [];
    for (let i = 0; i < 3; i++) for (let d = 0; d < 3; d++) T.push(L[dims[d].axis][i]);
    return { T, lattice: false };
  }
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

/** An axis in reciprocal lattice units (frame HKL, or r.l.u. without a frame). */
export const isHKL = (dim) => dim.frame === 'HKL' || (!dim.frame && /r\.?l\.?u/i.test(dim.units));

/** A real-space axis along a direct lattice vector (axis 0, 1 or 2 for a, b or c), in Å. */
export const isDirect = (dim) => dim.frame === 'direct' && [0, 1, 2].includes(dim.axis);

/**
 * Display geometry of the plane spanned by dims x and y: axis lengths per unit
 * coordinate, the cosine of the angle between them, and whether both axes
 * share a unit (equal aspect). Falls back to rectangular axes without a cell.
 */
export function planeGeometry(dims, cell, x, y) {
  const X = dims[x], Y = dims[y];
  if (cell && isDirect(X) && isDirect(Y)) {
    // The angle between two of a, b and c is the cell angle of the third.
    const angle = [cell.alpha, cell.beta, cell.gamma][3 - X.axis - Y.axis];
    return { lx: 1, ly: 1, cos: Math.cos(angle * Math.PI / 180), equal: true, lattice: true };
  }
  if (cell && X.basis && Y.basis && isHKL(X) && isHKL(Y)) {
    const G = reciprocalMetric(cell);
    const dot = (u, v) => u.reduce((s, ui, i) => s + ui * v.reduce((t, vj, j) => t + G[i][j] * vj, 0), 0);
    const lx = Math.sqrt(dot(X.basis.vec, X.basis.vec)), ly = Math.sqrt(dot(Y.basis.vec, Y.basis.vec));
    return { lx, ly, cos: dot(X.basis.vec, Y.basis.vec) / (lx * ly), equal: true, lattice: true };
  }
  return { lx: 1, ly: 1, cos: 0, equal: !!X.units && X.units === Y.units, lattice: false };
}
