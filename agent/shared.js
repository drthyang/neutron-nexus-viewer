// Fields and helpers the tools share: the volume id, outputs, axis lookup and value statistics.

import { z } from 'zod';

import { isDirect, isHKL } from '../js/nexus.js';
import { ToolError } from './tool.js';

export const volumeId = z.string().min(1).describe('The id open_volume returned.');
export const outputPath = z.string().min(1).describe('File to write; `~` is the home directory, relative paths start from the server\'s working directory.');
export const overwrite = z.boolean().optional().describe('Replace output_path if it exists (default false: an existing file is an error).');
export const vec3 = z.array(z.number()).length(3);

export const WRITTEN_BY = 'NeXus Viewer agent tools (https://github.com/drthyang/neutron-nexus-viewer)';

/** The index of a display axis given as its label ("L", "[H,H,0]", case-insensitive for letters) or index 0-2. */
export function axisIndex(dims, axis) {
  if (typeof axis === 'number') {
    if (![0, 1, 2].includes(axis)) throw new ToolError(`Axis index ${axis} is not 0, 1 or 2.`);
    return axis;
  }
  const a = String(axis).trim();
  let i = dims.findIndex((d) => d.label === a || d.name === a || d.longName === a);
  if (i < 0) i = dims.findIndex((d) => d.label.toLowerCase() === a.toLowerCase());
  if (i < 0) throw new ToolError(`No axis '${a}'. Axes: ${dims.map((d, k) => `${k} = ${d.label}`).join(', ')}.`, { axes: dims.map((d) => d.label) });
  return i;
}

/** Bin width of a display axis with uniform bins (the mean width otherwise). */
export const binWidth = (dim) => (dim.edges[dim.edges.length - 1] - dim.edges[0]) / (dim.edges.length - 1);

export const isUniform = (dim) => {
  const n = dim.edges.length - 1, w = binWidth(dim);
  return dim.edges.every((e, k) => Math.abs(e - dim.edges[0] - k * w) <= 1e-4 * w);
};

/** What kind of space the axes span. */
export function spaceOf(dims) {
  if (dims.every(isDirect)) return 'real space (Å along a, b, c: a 3D-ΔPDF)';
  if (dims.every(isHKL)) return 'reciprocal space (r.l.u.)';
  return 'other axes';
}

/** Name of the plane normal to display axis `fixed`, from the other two labels: "HK", or "[H,H,0]–L". */
export function planeName(dims, fixed) {
  const [x, y] = [0, 1, 2].filter((d) => d !== fixed).map((d) => dims[d].label);
  return x.length === 1 && y.length === 1 ? `${x}${y}` : `${x}–${y}`;
}

/** Quantiles of the finite values of an array, from at most `samples` evenly spread entries. */
export function valueStats(values, { samples = 2_000_000, positiveOnly = false } = {}) {
  const stride = Math.max(1, Math.floor(values.length / samples));
  const picked = [];
  let sum = 0, count = 0, min = Infinity, max = -Infinity;
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (v !== v) continue;
    count++;
    sum += v;
    if (v < min) min = v;
    if (v > max) max = v;
    if (i % stride === 0 && (!positiveOnly || v > 0)) picked.push(v);
  }
  if (!count) return { count: 0 };
  const sorted = Float64Array.from(picked).sort();
  const q = (p) => sorted[Math.min(sorted.length - 1, Math.floor(p * sorted.length))];
  return { count, min, max, mean: sum / count, median: q(0.5), p97: q(0.97), ...(stride > 1 ? { quantiles_from: `every ${stride}th voxel` } : {}) };
}

/** Rounds to `n` significant figures (finite values only). */
export const sig = (x, n = 6) => (Number.isFinite(x) && x !== 0 ? Number(x.toPrecision(n)) : x);
