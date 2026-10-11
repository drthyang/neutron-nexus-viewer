// The volumes a session has opened, with the processing the viewer would apply to them: a symmetry
// group (pooled in slices, cuts, I(Q), masks and the export) and a mask (coverage-edge erosion and
// symmetry outliers). Each volume stays in memory until it is closed, as in a viewer tab; the tools
// call the viewer's own modules (js/), so for the same choices they give the viewer's numbers.

import { closeSync, openSync, readSync, statSync } from 'node:fs';
import { basename } from 'node:path';

import h5wasm from 'h5wasm/node';

import { describeFile, isDirect, loadVariance, loadVolume } from '../js/nexus.js';
import { IDENTITY_MAP } from '../js/slab.js';
import { closeGroup, formatOp, indexMaps, metricChange, parseOps, PRESETS } from '../js/symmetry.js';
import { resolvePath } from './files.js';
import { ToolError } from './tool.js';

const SIGNATURE = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a];
const NO_SYMMETRY = { name: '1', ops: [[1, 0, 0, 0, 1, 0, 0, 0, 1]], maps: [IDENTITY_MAP] };

/** Whether `file` starts with the HDF5 signature (after a user block of 512, 1024 or 2048 bytes). */
function isHDF5(file) {
  const head = Buffer.alloc(2048 + 8), fd = openSync(file, 'r');
  try {
    readSync(fd, head, 0, head.length, 0);
  } finally {
    closeSync(fd);
  }
  return [0, 512, 1024, 2048].some((at) => SIGNATURE.every((b, i) => head[at + i] === b));
}

const groupKey = (ops) => ops.map((m) => m.join()).sort().join(';');
const PRESET_KEYS = new Map(PRESETS.map(([name, gens]) => [groupKey(closeGroup(parseOps(gens))), name]));

/** One open volume: the file's description (info), the float32 volume, and its processing. */
export class Dataset {
  constructor(id, path, info, loaded, seconds) {
    Object.assign(this, { id, path, name: basename(path), info, volume: loaded.volume, stats: loaded.stats, seconds });
    this.variance = undefined; // read when first needed; null without uncertainties
    this.varianceNote = '';
    this.symmetry = NO_SYMMETRY;
    this.mask = null; // { bits, radius, k, group, stats }
  }

  /** σ² per voxel, read from the file the first time; null when the file has none. */
  readVariance() {
    if (this.variance !== undefined) return this.variance;
    this.variance = null;
    if (!this.info.errors) return null;
    const file = new h5wasm.File(this.path, 'r');
    try {
      this.variance = loadVariance(file, this.info);
    } catch (err) {
      this.varianceNote = `uncertainties not read: ${err.message}`;
    } finally {
      file.close();
    }
    return this.variance;
  }

  /**
   * Set the symmetry group from operations ("k,-h+k,l; -h,-k,-l"; generators are closed into a group)
   * or a Laue class preset name. Returns what was set, with the metric check the viewer shows.
   */
  setSymmetry({ preset, operations }) {
    let text = operations;
    if (preset !== undefined) {
      const hit = PRESETS.find(([name]) => name === preset.trim());
      if (!hit) throw new ToolError(`No Laue class preset "${preset}".`, { presets: PRESETS.map(([name]) => name) });
      text = hit[1];
    }
    let ops, maps;
    try {
      ops = closeGroup(parseOps(text ?? ''));
      maps = indexMaps(ops, this.info.dims);
    } catch (err) {
      throw new ToolError(err.message);
    }
    const name = PRESET_KEYS.get(groupKey(ops)) ?? 'custom';
    this.symmetry = ops.length === 1 ? NO_SYMMETRY : { name, ops, maps };
    return this.describeSymmetry();
  }

  /** The current group: its name, operations and how well they keep the cell's metric. */
  describeSymmetry() {
    const { name, ops } = this.symmetry, { dims, lattice } = this.info;
    const out = { group: name, order: ops.length, operations: ops.map(formatOp).join('; ') };
    if (ops.length > 1) {
      const hkl = dims.every((d) => d.basis) || dims.every((d) => isDirect(d) && d.length > 0);
      if (lattice && hkl) {
        const change = metricChange(ops, lattice);
        out.metric_change = change;
        if (change > 0.02) out.warning = `These operations change the cell metric by up to ${(100 * change).toFixed(0)} %: they are not symmetries of this lattice. Check the setting.`;
      } else if (!hkl) {
        out.note = 'No HKL basis: the operations act on the display axes as they are.';
      }
    }
    return out;
  }

  /** The current group's index maps on another grid (the export's padded one). */
  mapsOn(dims) {
    if (this.symmetry.ops.length === 1) return [IDENTITY_MAP];
    try {
      return indexMaps(this.symmetry.ops, dims);
    } catch (err) {
      throw new ToolError(err.message);
    }
  }

  /** The processing in effect, for every result. */
  processing() {
    return {
      symmetry: this.symmetry.ops.length > 1 ? `${this.symmetry.name} (${this.symmetry.ops.length} operations)` : 'none',
      mask: this.mask ? maskSummary(this.mask) : 'none',
    };
  }
}

export function maskSummary({ radius, k, group, stats }) {
  const parts = [];
  if (radius) parts.push(`coverage-edge erosion ${radius} voxel${radius === 1 ? '' : 's'}`);
  if (k) parts.push(`outlier cut ${k}σ among ${group} equivalents`);
  const removed = stats.measured ? (stats.edge + stats.outlier) / stats.measured : 0;
  return `${parts.join(', ')}: ${(100 * removed).toFixed(2)} % of measured voxels removed`;
}

export class Workspace {
  #volumes = new Map();

  /** Open the 3-D histogram in `path` as `id` (default: the file name); `alignTo` displays its axes in another volume's order. */
  async open(path, { id, alignTo } = {}) {
    const file = resolvePath(path);
    let size;
    try {
      size = statSync(file).size;
    } catch {
      throw new ToolError(`No file at ${file}.`);
    }
    if (!size || !isHDF5(file)) throw new ToolError(`${file} is not an HDF5/NeXus file.`);
    const order = alignTo !== undefined ? this.get(alignTo).info.dims.map((d) => d.label) : null;
    await h5wasm.ready;
    const t0 = performance.now();
    let h5;
    try {
      h5 = new h5wasm.File(file, 'r');
    } catch {
      throw new ToolError(`${file} could not be opened as an HDF5/NeXus file.`);
    }
    let info, loaded;
    try {
      info = describeFile(h5, { order });
      loaded = loadVolume(h5, info);
    } catch (err) {
      throw new ToolError(err.message);
    } finally {
      h5.close();
    }
    const name = id ?? this.#freeId(basename(file).replace(/\.(nxs|h5|hdf5?|nx5)$/i, '').replace(/[^\w.-]+/g, '_'));
    const dataset = new Dataset(name, file, info, loaded, (performance.now() - t0) / 1000);
    this.#volumes.set(name, dataset);
    return dataset;
  }

  #freeId(stem) {
    let id = stem || 'volume';
    for (let k = 2; this.#volumes.has(id); k++) id = `${stem}_${k}`;
    return id;
  }

  /** The open volume `id`, or an error listing the open ones. */
  get(id) {
    const d = this.#volumes.get(id);
    if (!d) {
      const open = [...this.#volumes.keys()];
      throw new ToolError(`No volume '${id}' is open. ${open.length ? `Open: ${open.join(', ')}.` : 'Open one with open_volume.'}`, { open });
    }
    return d;
  }

  close(id) {
    this.get(id);
    this.#volumes.delete(id);
  }

  list() {
    return [...this.#volumes.values()];
  }
}
