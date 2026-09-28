// Module worker: opens the HDF5 file with h5wasm (random access through
// WORKERFS, so the file is never copied into memory), loads the histogram
// into a float32 volume and answers slice requests from the page.

import h5wasm from 'https://cdn.jsdelivr.net/npm/h5wasm@0.10.3/dist/esm/hdf5_hl.js';
import { binVolume, coarseGrid, orbitMean, surfaceNets } from './iso.js';
import { describeFile, loadVolume } from './nexus.js';
import { averageSlab, selectBins } from './slab.js';
import { indexMaps } from './symmetry.js';

let info = null, volume = null;
const coarse = new Map(), means = new Map();

self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'open') await open(data.file);
    else if (data.type === 'slice') slice(data);
    else if (data.type === 'iso') iso(data);
  } catch (err) {
    self.postMessage({ type: data.type === 'iso' ? 'iso-error' : 'error', id: data.id, fixed: data.fixed, message: err?.message ?? String(err) });
  }
};

// HDF5 superblock signature; it may follow a user block of 512, 1024, 2048... bytes.
const SIGNATURE = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a];

async function isHDF5(file) {
  const head = new Uint8Array(await file.slice(0, 2048 + 8).arrayBuffer());
  return [0, 512, 1024, 2048].some((at) => SIGNATURE.every((b, i) => head[at + i] === b));
}

async function open(file) {
  self.postMessage({ type: 'progress', label: 'Starting HDF5 reader', fraction: 0 });
  if (!(await isHDF5(file))) throw new Error(`${file.name} is not an HDF5/NeXus file.`);
  const { FS } = await h5wasm.ready;
  FS.mkdir('/work');
  FS.mount(FS.filesystems.WORKERFS, { files: [file] }, '/work');
  let h5;
  try {
    h5 = new h5wasm.File(`/work/${file.name}`, 'r');
  } catch {
    throw new Error(`${file.name} could not be opened as an HDF5/NeXus file.`);
  }
  try {
    info = describeFile(h5);
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
    coarse.set(maxBins, { grid, binned: binVolume(volume, info.shape, grid) });
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
function slice({ id, fixed, center, thickness, maps, symmetry }) {
  if (!volume) throw new Error('No histogram loaded.');
  const edges = info.dims[fixed].edges;
  const ids = selectBins(edges, center, thickness);
  if (!ids.length) throw new Error('No bins selected: increase the thickness or move the center.');
  const result = averageSlab(volume, info.shape, fixed, ids, maps);
  self.postMessage({
    type: 'slice', id, fixed, center, thickness, symmetry, order: maps.length, ...result,
    bins: ids.length, slab: [edges[ids[0]], edges[ids.at(-1) + 1]],
  }, [result.values.buffer, result.counts.buffer]);
}
