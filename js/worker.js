// Module worker: opens the HDF5 file with h5wasm (random access through
// WORKERFS, so the file is never copied into memory), loads the histogram
// into a float32 volume and answers slice requests from the page.

import h5wasm from 'https://cdn.jsdelivr.net/npm/h5wasm@0.10.3/dist/esm/hdf5_hl.js';
import { lineCut } from './cut.js';
import { symmetrizeForExport, writeNebulaFile } from './export.js';
import { binVolume, coarseGrid, orbitMean, surfaceNets } from './iso.js';
import { edgeMask, maskStats, outlierMask } from './mask.js';
import { describeFile, loadVariance, loadVolume } from './nexus.js';
import { powderAverage } from './powder.js';
import { measuredBox, reindexGeometry, reindexVolume } from './reindex.js';
import { averageSlab, selectBins } from './slab.js';
import { indexMaps } from './symmetry.js';

// `variance` (σ² per voxel) is read from the file when I(Q) or a line cut first needs it:
// undefined until then, null without uncertainties.
let info = null, volume = null, mask = null, path = null, variance, varianceNote = '';
const coarse = new Map(), means = new Map();

self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'open') await open(data.file, data.order);
    else if (data.type === 'slice') slice(data);
    else if (data.type === 'iso') iso(data);
    else if (data.type === 'mask') buildMask(data);
    else if (data.type === 'mask-download') await downloadMask(data);
    else if (data.type === 'export') await exportVolume(data);
    else if (data.type === 'powder') powder(data);
    else if (data.type === 'cut') cut(data);
    else if (data.type === 'reindex-plan') reindexPlan(data);
    else if (data.type === 'reindex') await reindex(data);
  } catch (err) {
    const type = {
      iso: 'iso-error', mask: 'mask-error', 'mask-download': 'mask-error', export: 'export-error', powder: 'powder-error', cut: 'cut-error',
      'reindex-plan': 'reindex-error', reindex: 'reindex-error',
    }[data.type] ?? 'error';
    self.postMessage({ type, id: data.id, fixed: data.fixed, message: err?.message ?? String(err) });
  }
};

// HDF5 superblock signature; it may follow a user block of 512, 1024, 2048... bytes.
const SIGNATURE = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a];

async function isHDF5(file) {
  const head = new Uint8Array(await file.slice(0, 2048 + 8).arrayBuffer());
  return [0, 512, 1024, 2048].some((at) => SIGNATURE.every((b, i) => head[at + i] === b));
}

// `order`: the labels of another dataset's display axes, to display the same axes in that order.
async function open(file, order = null) {
  self.postMessage({ type: 'progress', label: 'Starting HDF5 reader', fraction: 0 });
  if (!(await isHDF5(file))) throw new Error(`${file.name} is not an HDF5/NeXus file.`);
  const { FS } = await h5wasm.ready;
  FS.mkdir('/work');
  FS.mount(FS.filesystems.WORKERFS, { files: [file] }, '/work');
  path = `/work/${file.name}`;
  let h5;
  try {
    h5 = new h5wasm.File(path, 'r');
  } catch {
    throw new Error(`${file.name} could not be opened as an HDF5/NeXus file.`);
  }
  try {
    info = describeFile(h5, { order });
    self.postMessage({ type: 'meta', info });
    const t0 = performance.now();
    let last = 0;
    const loaded = loadVolume(h5, info, (fraction) => {
      const now = performance.now();
      if (now - last > 100 || fraction === 1) {
        last = now;
        self.postMessage({ type: 'progress', label: 'Reading histogram', fraction });
      }
    });
    volume = loaded.volume;
    self.postMessage({ type: 'ready', stats: loaded.stats, seconds: (performance.now() - t0) / 1000 });
  } finally {
    h5.close();
  }
}

// Isosurface of the binned (and symmetrized) volume. `ops` are the group's
// hkl matrices; their maps are rebuilt for the coarse grid. Without a level,
// the 99.5th percentile of positive block means is used.
function iso({ id, maxBins, level, ops, symmetry }) {
  if (!volume) throw new Error('No histogram loaded.');
  if (!coarse.has(maxBins)) {
    const grid = coarseGrid(info.dims, maxBins);
    coarse.set(maxBins, { grid, binned: binVolume(volume, info.shape, grid, mask) });
  }
  const { grid, binned } = coarse.get(maxBins);
  const key = `${maxBins}|${ops.map((m) => m.join()).join(';')}`;
  if (!means.has(key)) {
    let maps, note = '';
    try {
      maps = indexMaps(ops, grid.dims);
    } catch (err) {
      note = `3-D view not symmetrized: ${err.message}`;
      maps = indexMaps([ops[0]], grid.dims);
    }
    const mean = orbitMean(binned, grid.shape, maps);
    const positive = Float32Array.from(mean.filter((v) => v > 0)).sort();
    const quantile = (q) => positive[Math.floor(q * (positive.length - 1))] ?? 1;
    means.clear();
    means.set(key, { mean, note, range: [quantile(0.5), quantile(0.9999)], auto: quantile(0.995) });
  }
  const { mean, range, auto, note } = means.get(key);
  const value = level ?? auto;
  const { positions, indices } = surfaceNets(mean, grid.shape, value);
  // Grid indices -> display coordinates at block centers.
  for (let v = 0; v < positions.length; v += 3) {
    for (let d = 0; d < 3; d++) {
      const e = grid.dims[d].edges, w = e[1] - e[0];
      positions[v + d] = e[0] + (positions[v + d] + 0.5) * w;
    }
  }
  self.postMessage({
    type: 'iso', id, level: value, range, auto, note, symmetry, factor: grid.factor, shape: grid.shape,
    positions, indices,
  }, [positions.buffer, indices.buffer]);
}

// `maps` are the index maps of the full symmetry group (identity included);
// `symmetry` is its display name, echoed back with the result.
// With `removed`, only voxels removed by the user mask are averaged.
function slice({ id, fixed, center, thickness, maps, symmetry, removed }) {
  if (!volume) throw new Error('No histogram loaded.');
  const edges = info.dims[fixed].edges;
  const ids = selectBins(edges, center, thickness);
  if (!ids.length) throw new Error('No bins selected: increase the thickness or move the center.');
  const result = averageSlab(volume, info.shape, fixed, ids, maps, mask, !!(removed && mask));
  self.postMessage({
    type: 'slice', id, fixed, center, thickness, symmetry, order: maps.length, removed: !!(removed && mask), ...result,
    bins: ids.length, slab: [edges[ids[0]], edges[ids.at(-1) + 1]],
  }, [result.values.buffer, result.counts.buffer]);
}

// User mask on the unsymmetrized volume: coverage-edge erosion by `radius`
// voxels, then rejection of voxels more than k robust sigmas above the median
// of their symmetry equivalents (`maps`, fine grid). Zero disables either.
function buildMask({ id, radius, k, maps, symmetry }) {
  if (!volume) throw new Error('No histogram loaded.');
  const t0 = performance.now();
  let next = null;
  if (radius > 0 || k > 0) {
    self.postMessage({ type: 'progress-mask', label: 'Eroding coverage edges', fraction: 0 });
    next = edgeMask(volume, info.shape, radius);
    if (k > 0) {
      if (maps.length < 3) throw new Error('Outlier rejection needs a symmetry group with at least 3 operations.');
      outlierMask(volume, info.shape, maps, k, next, 3,
        (fraction) => self.postMessage({ type: 'progress-mask', label: `Comparing ${symmetry} equivalents`, fraction }));
    }
  }
  mask = next;
  coarse.clear();
  means.clear();
  const stats = mask ? maskStats(volume, mask) : null;
  self.postMessage({ type: 'mask', id, stats, radius, k, symmetry, seconds: (performance.now() - t0) / 1000 });
}

// The mask as a gzipped NumPy .npy (uint8, storage order like the signal
// dataset; 1 = coverage edge, 2 = symmetry outlier).
async function downloadMask() {
  if (!mask) throw new Error('No mask to download.');
  const dict = `{'descr': '|u1', 'fortran_order': False, 'shape': (${info.shape.join(', ')}), }`;
  const padded = dict + ' '.repeat(63 - ((10 + dict.length) % 64)) + '\n';
  const header = new Uint8Array(10 + padded.length);
  header.set([0x93, ...new TextEncoder().encode('NUMPY'), 1, 0, padded.length & 255, padded.length >> 8]);
  header.set(new TextEncoder().encode(padded), 10);
  const stream = new Blob([header, mask]).stream().pipeThrough(new CompressionStream('gzip'));
  self.postMessage({ type: 'mask-file', blob: await new Response(stream).blob() });
}

// I(Q), the powder average of the masked volume (see powder.js), on the shells
// of `plan` (powderPlan() on the page) with the group's `maps` on this grid.
function powder({ id, plan, maps, symmetry }) {
  if (!volume) throw new Error('No histogram loaded.');
  const t0 = performance.now();
  let last = 0;
  const progress = (label, from, to) => (fraction) => {
    const now = performance.now();
    if (now - last > 100 || fraction === 1) {
      last = now;
      self.postMessage({ type: 'progress-powder', label, fraction: from + (to - from) * fraction });
    }
  };
  // The first time, a fifth of the progress is reading the uncertainties.
  const start = variance === undefined && info.errors ? 0.2 : 0;
  readVariance(progress('Reading uncertainties', 0, start));
  const result = powderAverage(volume, info.shape, plan, maps, mask, variance, progress('Averaging shells', start, 1));
  const errors = variance ? info.errors.path.split('/').pop() : null;
  self.postMessage({
    type: 'powder', id, ...result, bins: plan.bins, split: plan.split, frame: plan.frame, symmetry, order: maps.length,
    masked: !!mask, errors, errorsNote: varianceNote, seconds: (performance.now() - t0) / 1000,
  }, [result.edges.buffer, result.intensity.buffer, result.sigma.buffer, result.voxels.buffer, result.coverage.buffer]);
}

/** Read the uncertainties the first time they are needed; `variance` stays null without them. */
function readVariance(onProgress) {
  if (variance !== undefined) return;
  variance = null;
  if (!info.errors) return;
  const h5 = new h5wasm.File(path, 'r');
  try {
    variance = loadVariance(h5, info, onProgress);
  } catch (err) {
    varianceNote = `uncertainties not read: ${err.message}`;
  } finally {
    h5.close();
  }
}

// A line cut through the slab of a slice (see cut.js); `cut` is its geometry in
// display coordinates, the same for both datasets, and `maps` the group on this grid.
function cut({ id, cut: spec, maps, symmetry, removed }) {
  if (!volume) throw new Error('No histogram loaded.');
  const t0 = performance.now();
  let last = 0;
  readVariance((fraction) => {
    const now = performance.now();
    if (now - last > 100 || fraction === 1) {
      last = now;
      self.postMessage({ type: 'progress-cut', label: 'Reading uncertainties', fraction });
    }
  });
  const invert = !!(removed && mask);
  const result = lineCut(volume, info.shape, info.dims, spec, maps, mask, invert, variance);
  const errors = variance ? info.errors.path.split('/').pop() : null;
  self.postMessage({
    type: 'cut', id, ...result, symmetry, order: maps.length, masked: !!mask, removed: invert, errors, errorsNote: varianceNote,
    seconds: (performance.now() - t0) / 1000,
  }, [result.edges.buffer, result.intensity.buffer, result.sigma.buffer, result.voxels.buffer]);
}

// Reindexing with another UB (see reindex.js): `T` maps the file's indices to the new ones.
// The plan is the box of the measured voxels in the new indices.
function reindexPlan({ id, T }) {
  if (!volume) throw new Error('No histogram loaded.');
  const { TW, axes } = reindexGeometry(info.dims, T);
  const box = measuredBox(volume, info.shape, axes, TW);
  if (!box) throw new Error('The volume has no measured voxels.');
  self.postMessage({ type: 'reindex-plan', id, box });
}

// The volume as measured (the user mask and symmetry are not applied), resampled onto
// `grid` with nsub³ sub-samples per voxel, as a Mantid MDHistoWorkspace with the new UB.
async function reindex({ id, T, grid, nsub, ub, cell, title, logs }) {
  if (!volume) throw new Error('No histogram loaded.');
  const t0 = performance.now();
  let last = 0;
  const progress = (label, from, to) => (fraction) => {
    const now = performance.now();
    if (now - last > 100 || fraction === 1) {
      last = now;
      self.postMessage({ type: 'progress-reindex', id, label, fraction: from + (to - from) * fraction });
    }
  };
  const start = variance === undefined && info.errors ? 0.1 : 0;
  readVariance(progress('Reading uncertainties', 0, start));
  const { M, axes } = reindexGeometry(info.dims, T);
  const out = reindexVolume(volume, variance, info.shape, axes, M, grid, nsub, progress('Resampling', start, 0.8));
  self.postMessage({ type: 'progress-reindex', id, label: 'Writing HDF5', fraction: 0.8 });
  const [{ FS }, { writeMantidMD }] = await Promise.all([h5wasm.ready, import('./rigaku-reduce.js')]);
  const file = '/reindexed.nxs';
  writeMantidMD(h5wasm, file, {
    shape: [grid.shape[2], grid.shape[1], grid.shape[0]], edges: grid.edges, signal: out.signal, errors2: out.errors2, events: out.events,
    ub, cell, title, logs,
  });
  const blob = new Blob([FS.readFile(file)], { type: 'application/x-hdf5' });
  FS.unlink(file);
  self.postMessage({
    type: 'reindex-file', id, blob, covered: out.covered, voxels: grid.voxels, errors: !!variance, errorsNote: varianceNote,
    seconds: (performance.now() - t0) / 1000,
  });
}

// The symmetrized, masked volume as a NEBULA3D input file (see export.js):
// `plan` from exportPlan() on the page, `maps` on its padded grid.
async function exportVolume({ id, plan, maps, attrs }) {
  if (!volume) throw new Error('No histogram loaded.');
  const t0 = performance.now();
  let last = 0;
  const result = symmetrizeForExport(volume, info.shape, plan, maps, mask, (fraction) => {
    const now = performance.now();
    if (now - last > 100 || fraction === 1) {
      last = now;
      self.postMessage({ type: 'progress-export', label: 'Symmetrizing', fraction: 0.75 * fraction });
    }
  });
  self.postMessage({ type: 'progress-export', label: 'Writing HDF5', fraction: 0.75 });
  const { FS } = await h5wasm.ready;
  const path = '/nebula3d-export.nxs';
  writeNebulaFile(h5wasm, path, plan, result, attrs);
  const blob = new Blob([FS.readFile(path)], { type: 'application/x-hdf5' });
  FS.unlink(path);
  self.postMessage({ type: 'export-file', id, blob, stats: result.stats, seconds: (performance.now() - t0) / 1000 });
}
