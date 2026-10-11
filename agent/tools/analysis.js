// What the viewer shows, as numbers (and images): slices, line cuts, I(Q) and the value at points,
// each with the volume's symmetry averaging and mask, computed by the viewer's own modules.

import { readFileSync } from 'node:fs';

import { z } from 'zod';

import { cutAxis, cutEdges, lineCut } from '../../js/cut.js';
import { cartesianBasis, nominalCell, planeGeometry } from '../../js/nexus.js';
import { parseBins, powderAverage, powderPlan } from '../../js/powder.js';
import { averageSlab, selectBins } from '../../js/slab.js';
import { npyBytes, num, resolvePath, writeOutput } from '../files.js';
import { COLORMAP_NAMES, renderSlice } from '../png.js';
import { axisIndex, binWidth, outputPath, overwrite, planeName, sig, valueStats, vec3, volumeId, WRITTEN_BY } from '../shared.js';
import { defineTool, ToolError, withImages } from '../tool.js';

const anglesField = z.enum(['nominal', 'measured']).optional()
  .describe('Cell angles for distances and drawing: nominal (default, as the viewer: direct-cell angles within 1° of 60/90/120° snapped) or measured.');
const removedOnly = z.boolean().optional().describe('Average only the voxels the mask removed (to inspect what it takes out).');

/** The cell for drawing and distances, as the viewer's "Nominal cell angles" switch. */
const displayCell = (lattice, angles) => lattice && (angles === 'measured' ? lattice : nominalCell(lattice));

const maskArgs = (d, removed) => {
  if (removed && !d.mask) throw new ToolError('removed_only needs a mask: build one with build_mask.');
  return [d.mask?.bits ?? null, !!removed];
};

const header = (d, lines) => [
  `# Written by ${WRITTEN_BY} on ${new Date().toISOString()}`,
  `# Volume: ${d.path}`,
  `# Symmetry: ${d.processing().symmetry}${d.symmetry.ops.length > 1 ? ', equivalent voxels pooled' : ''}`,
  `# Mask: ${d.processing().mask}`,
  ...lines,
];

// ---------------------------------------------------------------- slices

/** The viewer's Auto color range for one slice: vmin 0 (or -vmax for signed data), vmax at the 97th percentile, asinh softening at the median. */
function autoColor(values, signed, scale) {
  const positive = [];
  for (const v of values) if (v > 0 || (signed && v < 0)) positive.push(Math.abs(v));
  positive.sort((a, b) => a - b);
  const quantile = (q) => positive[Math.min(positive.length - 1, Math.floor(q * positive.length))];
  const vmax = positive.length ? sig(quantile(0.97), 2) : 1;
  return {
    min: scale === 'log' ? sig(vmax / 1000, 2) : signed ? -vmax : 0,
    max: vmax,
    soft: positive.length ? sig(quantile(0.5), 1) : sig(vmax / 20, 2),
  };
}

export const sliceVolume = defineTool({
  name: 'slice_volume',
  title: 'Slice a volume',
  description:
    'A 2-D slice of an open volume, as the viewer draws it: the slab of bins whose centres lie within thickness/2 of center along fixed_axis, averaged with the volume\'s symmetry (each pixel pools the distinct measured, unmasked voxels of the orbits in its slab column) and mask. ' +
    'Returns the slab, axes, value statistics and the viewer\'s Auto color range, and by default a PNG of the slice in its true geometry (oblique axes at their real angle, no smoothing; grey = no data) to look at. ' +
    'output_path writes the values as .npy (float32, rows = y bins, columns = x bins, NaN without data), .csv (x, y, value, voxels per pixel) or .png.',
  input: z.strictObject({
    volume_id: volumeId,
    fixed_axis: z.union([z.string(), z.number().int()]).describe('The axis held fixed, by label (e.g. "L" for the HK plane) or index 0-2.'),
    center: z.number().optional().describe('Slab centre along fixed_axis (default 0 when in range, else the middle).'),
    thickness: z.number().positive().optional().describe('Slab thickness along fixed_axis (default one bin).'),
    removed_only: removedOnly,
    image: z.boolean().optional().describe('Return the PNG (default true).'),
    color: z.strictObject({
      scale: z.enum(['asinh', 'linear', 'log']).optional().describe('Default asinh, as the viewer.'),
      min: z.number().optional(),
      max: z.number().optional(),
      soft: z.number().positive().optional().describe('Asinh softening: nearly linear below it, logarithmic above.'),
      colormap: z.enum(COLORMAP_NAMES).optional().describe('Default magma (coolwarm for a ΔPDF).'),
    }).optional().describe('Color scale of the image; unset limits are the viewer\'s Auto ones.'),
    angles: anglesField,
    image_size: z.number().int().min(64).max(1600).optional().describe('Longest side of the image in pixels (default 640).'),
    output_path: outputPath.optional(),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  run(args, { workspace }) {
    const d = workspace.get(args.volume_id), { dims, shape, lattice, signed } = d.info;
    const fixed = axisIndex(dims, args.fixed_axis), dim = dims[fixed];
    const [x, y] = [0, 1, 2].filter((k) => k !== fixed);
    const lo = dim.edges[0], hi = dim.edges[dim.edges.length - 1];
    const center = args.center ?? (lo <= 0 && hi >= 0 ? 0 : (lo + hi) / 2);
    const thickness = args.thickness ?? binWidth(dim);
    const ids = selectBins(dim.edges, center, thickness);
    if (!ids.length) throw new ToolError(`No bins of ${dim.label} lie within ${thickness / 2} of ${center} (range ${lo} to ${hi}): increase the thickness or move the center.`);
    const [mask, removed] = maskArgs(d, args.removed_only);
    const r = averageSlab(d.volume, shape, fixed, ids, d.symmetry.maps, mask, removed);
    const scale = args.color?.scale ?? 'asinh';
    const auto = autoColor(r.values, signed, scale);
    const color = {
      scale,
      min: args.color?.min ?? auto.min,
      max: args.color?.max ?? auto.max,
      soft: args.color?.soft ?? auto.soft,
      colormap: args.color?.colormap ?? (signed ? 'coolwarm' : 'magma'),
    };
    if (!(color.max > color.min)) throw new ToolError(`The color range needs max (${color.max}) above min (${color.min}).`);
    if (scale === 'log' && !(color.min > 0)) throw new ToolError('A log scale needs a positive min.');
    const geometry = planeGeometry(dims, displayCell(lattice, args.angles), x, y);
    const angle = Math.acos(Math.max(-1, Math.min(1, geometry.cos))) * 180 / Math.PI;
    const axis = (k) => ({ label: dims[k].label, min: dims[k].edges[0], max: dims[k].edges[dims[k].edges.length - 1], bins: dims[k].edges.length - 1 });
    const out = {
      volume_id: d.id,
      plane: planeName(dims, fixed),
      slab: { axis: dim.label, center, thickness, bins: ids.length, range: [dim.edges[ids[0]], dim.edges[ids[ids.length - 1] + 1]] },
      x_axis: axis(x),
      y_axis: axis(y),
      ...(geometry.lattice ? { angle_between_axes_deg: angle } : {}),
      shape: [r.rows, r.cols],
      coverage: r.coverage,
      values: valueStats(r.values),
      color,
      processing: { ...d.processing(), ...(removed ? { removed_only: true } : {}) },
    };
    const images = [];
    const wantPng = args.image !== false || /\.png$/i.test(args.output_path ?? '');
    if (wantPng) {
      const img = renderSlice({ values: r.values, rows: r.rows, cols: r.cols, xEdges: dims[x].edges, yEdges: dims[y].edges, geometry, color, colormap: color.colormap, maxSize: args.image_size ?? 640 });
      out.image = {
        width: img.width,
        height: img.height,
        orientation: `${dims[x].label} increases to the right, ${dims[y].label} upward${geometry.lattice ? ` at ${sig(angle, 4)}° to ${dims[x].label}` : ''}; grey = no data, white = outside the volume`,
      };
      if (args.image !== false) images.push({ data: img.png, mimeType: 'image/png' });
      if (/\.png$/i.test(args.output_path ?? '')) out.written = writeOutput(args.output_path, img.png, args.overwrite);
    }
    if (args.output_path && !/\.png$/i.test(args.output_path)) {
      if (/\.npy(\.gz)?$/i.test(args.output_path)) {
        out.written = writeOutput(args.output_path, npyBytes(r.values, [r.rows, r.cols]), args.overwrite);
      } else if (/\.(csv|txt|dat)$/i.test(args.output_path)) {
        const cx = centers(dims[x].edges), cy = centers(dims[y].edges);
        const lines = header(d, [
          `# ${out.plane} slice, ${dim.label} = ${center} ± ${thickness / 2} (${ids.length} bins, ${out.slab.range.join(' to ')})`,
          `# value: equal-weight mean of the distinct measured, unmasked voxels in the symmetry orbits of the slab column; voxels: how many`,
          `${dims[x].label},${dims[y].label},value,voxels`,
        ]);
        for (let row = 0; row < r.rows; row++) for (let col = 0; col < r.cols; col++) {
          const k = row * r.cols + col;
          lines.push(`${num(cx[col])},${num(cy[row])},${num(r.values[k])},${r.counts[k]}`);
        }
        out.written = writeOutput(args.output_path, `${lines.join('\n')}\n`, args.overwrite);
      } else {
        throw new ToolError('output_path must end in .npy, .csv, .txt, .dat or .png.');
      }
    }
    return withImages(out, images);
  },
});

const centers = (edges) => Float64Array.from({ length: edges.length - 1 }, (_, k) => (edges[k] + edges[k + 1]) / 2);

// ---------------------------------------------------------------- line cuts

export const lineCutTool = defineTool({
  name: 'line_cut',
  title: 'Line cut through a volume',
  description:
    'The 1-D profile of an open volume along a straight line between two points in any direction in 3-D (e.g. (H, H, H) or along L through (0.5, 0, L)), as the viewer\'s Cut mode computes it: binned along the axis the line changes most along, each point the equal-weight mean of the distinct measured, unmasked voxels (with their symmetry orbits) whose centres lie within a rod of diameter `width` around the line and project within step/2 of the point. ' +
    'σ is propagated from the file\'s uncertainties (independent voxels). Points are given in the volume\'s axis coordinates (as its axes are labelled, e.g. [H, K, L]). output_path writes .csv (one row per point) or .txt (the viewer\'s commented text format).',
  input: z.strictObject({
    volume_id: volumeId,
    start: vec3.describe('One end, in the volume\'s axis coordinates (axis 0, 1, 2).'),
    end: vec3.describe('The other end.'),
    width: z.number().positive().optional().describe('Rod diameter: in Å⁻¹ with 2π when the file has a cell and hkl axes, otherwise in the axes\' unit (default: three voxels across, as the viewer).'),
    step: z.number().positive().optional().describe('Spacing of the points along the dominant axis (default its bin width).'),
    removed_only: removedOnly,
    angles: anglesField,
    limit: z.number().int().min(1).max(20000).optional().describe('Most points returned (default 400; output_path gets all).'),
    output_path: outputPath.optional(),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  run(args, { workspace }) {
    const d = workspace.get(args.volume_id), { dims, shape } = d.info;
    let [a, b] = [args.start, args.end];
    const dom = cutAxis(a, b);
    if (!(Math.abs(b[dom] - a[dom]) > 0)) throw new ToolError('The cut has no length: start and end are the same point.');
    const reversed = b[dom] < a[dom];
    if (reversed) [a, b] = [b, a];
    const { T, lattice } = cartesianBasis(dims, displayCell(d.info.lattice, args.angles));
    const scale = lattice ? 2 * Math.PI : 1;
    const unit = lattice ? 'Å⁻¹' : dims.every((x) => x.units === dims[0].units) ? dims[0].units : '';
    const voxel = Math.max(...dims.map((x, i) => binWidth(x) * Math.hypot(T[i], T[3 + i], T[6 + i])));
    const width = args.width ?? sig(3 * voxel * scale, 2);
    const step = args.step ?? binWidth(dims[dom]);
    let edges;
    try {
      edges = cutEdges(a, b, dom, step);
    } catch (err) {
      throw new ToolError(err.message);
    }
    const spec = { a, b, dom, edges, radius: width / scale / 2, T };
    const [mask, removed] = maskArgs(d, args.removed_only);
    const variance = d.readVariance();
    const r = lineCut(d.volume, shape, dims, spec, d.symmetry.maps, mask, removed, variance);
    const at = (s) => a.map((ad, k) => ad + ((s - a[dom]) / (b[dom] - a[dom])) * (b[k] - ad));
    const points = Array.from(r.intensity, (I, k) => {
      const s = (r.edges[k] + r.edges[k + 1]) / 2;
      return { [dims[dom].label]: s, at: at(s), intensity: I, sigma: r.sigma[k], voxels: r.voxels[k] };
    });
    const measured = points.filter((p) => p.voxels > 0);
    const peak = measured.reduce((m, p) => (p.intensity > (m?.intensity ?? -Infinity) ? p : m), null);
    const limit = args.limit ?? 400;
    const out = {
      volume_id: d.id,
      from: a,
      to: b,
      ...(reversed ? { note: 'start and end were swapped: a cut runs toward larger values of its axis.' } : {}),
      axis: dims[dom].label,
      step,
      width: { diameter: width, unit: unit || 'axis units' },
      points_total: points.length,
      points_with_data: measured.length,
      ...(peak ? { highest: { [dims[dom].label]: peak[dims[dom].label], intensity: peak.intensity, sigma: peak.sigma } } : {}),
      uncertainties: variance ? 'σ propagated from the file, voxels independent' : `none${d.varianceNote ? ` (${d.varianceNote})` : ' in the file'}`,
      processing: { ...d.processing(), ...(removed ? { removed_only: true } : {}) },
      points: points.slice(0, limit),
      ...(points.length > limit ? { truncated: `${points.length - limit} more points; give output_path for all` } : {}),
    };
    if (args.output_path) {
      const labels = dims.map((x) => x.label);
      let text;
      if (/\.csv$/i.test(args.output_path)) {
        text = [`${dims[dom].label},${labels.map((l) => `at_${l}`).join(',')},intensity,sigma,voxels`,
          ...points.map((p) => [p[dims[dom].label], ...p.at, p.intensity, p.sigma, p.voxels].map(num).join(','))].join('\n');
      } else {
        text = header(d, [
          `# Line cut from (${a.map(num).join(', ')}) to (${b.map(num).join(', ')}), points every ${num(step)} in ${dims[dom].label}`,
          `# Each point averages the voxels whose centres lie in a rod of diameter ${num(width)} ${unit || 'units'} around the line and project within ±${num(step / 2)} ${dims[dom].label} of it`,
          `# sigma: ${out.uncertainties}; voxels: distinct voxels with data pooled at the point`,
          `# ${dims[dom].label} ${labels.join(' ')} I sigma voxels`,
          ...points.map((p) => [p[dims[dom].label], ...p.at, p.intensity, p.sigma, p.voxels].map(num).join(' ')),
        ]).join('\n');
      }
      out.written = writeOutput(args.output_path, `${text}\n`, args.overwrite);
    }
    return out;
  },
});

// ---------------------------------------------------------------- I(Q)

export const powderAverageTool = defineTool({
  name: 'powder_average',
  title: 'I(Q): the powder average of a volume',
  description:
    'Reduces an open volume to I(|Q|), as the viewer\'s I(Q) view: per |Q| shell, the intensity integrated over the part of the shell with data divided by that part\'s volume (unmeasured and masked voxels are left out, not counted as zero), each symmetry orbit weighted by its multiplicity, voxels split into split³ sub-cells binned by their own |Q|, σ propagated from the file. ' +
    'Also gives each shell\'s coverage (fraction of its volume with data). Needs hkl axes and a unit cell, or Q axes in Å⁻¹. |Q| includes 2π. output_path writes the viewer\'s text format (.txt/.dat) or .csv.',
  input: z.strictObject({
    volume_id: volumeId,
    bins: z.string().optional().describe('Shells as Mantid Rebin parameters: one step ("0.02"; negative for logarithmic, "-0.01" = 1 %) or ranges "0.5,0.02,3,0.05,10" (default: the shortest bin step in |Q|).'),
    q_min: z.number().min(0).optional().describe('Start of a single-step range (default 0).'),
    q_max: z.number().positive().optional().describe('End of a single-step range (default the farthest grid corner).'),
    split: z.number().int().min(1).max(4).optional().describe('Sub-cells per axis per voxel (default 2, as the viewer).'),
    limit: z.number().int().min(1).max(20000).optional().describe('Most shells returned (default 400; output_path gets all).'),
    output_path: outputPath.optional(),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  run(args, { workspace }) {
    const d = workspace.get(args.volume_id), { dims, shape, lattice } = d.info;
    let plan;
    try {
      plan = powderPlan(dims, lattice, { bins: parseBins(args.bins ?? ''), qmin: args.q_min ?? null, qmax: args.q_max ?? null, split: args.split ?? 2 });
    } catch (err) {
      throw new ToolError(err.message);
    }
    const t0 = performance.now();
    const variance = d.readVariance();
    const r = powderAverage(d.volume, shape, plan, d.symmetry.maps, d.mask?.bits ?? null, variance);
    let last = r.voxels.length - 1;
    while (last > 0 && !(r.voxels[last] > 0)) last--;
    const shells = Array.from({ length: last + 1 }, (_, k) => ({
      q: (r.edges[k] + r.edges[k + 1]) / 2, intensity: r.intensity[k], sigma: r.sigma[k], coverage: r.coverage[k], voxels: r.voxels[k],
    }));
    const limit = args.limit ?? 400, first = Math.max(0, shells.findIndex((s) => s.voxels > 0));
    const out = {
      volume_id: d.id,
      q_unit: 'Å⁻¹, with 2π',
      q_from: plan.frame === 'HKL' ? `the cell (${lattice.source})` : 'the Q axes',
      shells_total: shells.length,
      shell_edges: { first: r.edges[0], last: r.edges[last + 1], step: plan.bins ?? 'given edges' },
      split: plan.split,
      orbits: r.orbits,
      uncertainties: variance ? 'σ propagated from the file' : `none${d.varianceNote ? ` (${d.varianceNote})` : ' in the file'}`,
      processing: d.processing(),
      seconds: (performance.now() - t0) / 1000,
      ...(first > 0 ? { empty_below: `the first ${first} shells (|Q| < ${num(r.edges[first])}) have no data and are not listed` } : {}),
      shells: shells.slice(first, first + limit),
      ...(shells.length - first > limit ? { truncated: `${shells.length - first - limit} more shells; give output_path for all` } : {}),
    };
    if (args.output_path) {
      const text = /\.csv$/i.test(args.output_path)
        ? ['q,intensity,sigma,coverage,voxels', ...shells.map((s) => [s.q, s.intensity, s.sigma, s.coverage, s.voxels].map(num).join(','))].join('\n')
        : header(d, [
          `# I(Q), the powder average; |Q| in 1/Angstrom with 2*pi from ${out.q_from}`,
          `# Shells: ${plan.bins ? `Rebin ${plan.bins.join(', ')}` : 'given edges'}; voxels split into ${plan.split}^3 sub-cells; Q is the shell centre`,
          '# I: intensity integrated over the part of the shell with data / volume of that part; coverage: fraction of the shell with data',
          `# sigma: ${out.uncertainties}`,
          '# Q I sigma coverage voxels',
          ...shells.map((s) => [s.q, s.intensity, s.sigma, s.coverage, s.voxels].map(num).join(' ')),
        ]).join('\n');
      out.written = writeOutput(args.output_path, `${text}\n`, args.overwrite);
    }
    return out;
  },
});

// ---------------------------------------------------------------- values at points

/** Points from a NEXPLAN bragg_positions file (JSON, schema nexplan-bragg-nodes/1, or its CSV). */
function readNodes(path, dims) {
  const file = resolvePath(path);
  let text;
  try {
    text = readFileSync(file, 'utf8');
  } catch {
    throw new ToolError(`Cannot read ${file}.`);
  }
  if (/^\s*\{/.test(text)) {
    const json = JSON.parse(text);
    if (!String(json.schema ?? '').startsWith('nexplan-bragg-nodes/')) throw new ToolError(`${file} is not a NEXPLAN bragg_positions file (schema nexplan-bragg-nodes/1).`);
    checkAxes(json.grid?.axes, dims, file);
    return json.nodes.map((n) => ({ at: n.grid, hkl: n.hkl, class: n.class }));
  }
  const rows = text.split(/\r?\n/).filter((l) => l && !l.startsWith('#'));
  const head = rows.shift()?.split(',') ?? [];
  const col = (name) => head.indexOf(name);
  if ([col('x0'), col('x1'), col('x2')].some((c) => c < 0)) throw new ToolError(`${file} has no x0, x1, x2 columns (grid coordinates).`);
  return rows.map((line) => {
    const v = line.split(',');
    return { at: [0, 1, 2].map((i) => Number(v[col(`x${i}`)])), hkl: ['h', 'k', 'l'].map((c) => Number(v[col(c)])), class: v[col('class')] };
  });
}

function checkAxes(axes, dims, file) {
  if (!axes) return;
  const same = dims.every((d, i) => d.basis && d.basis.vec.every((v, j) => Math.abs(v - axes[i][j]) < 1e-9));
  if (!same) throw new ToolError(`${file} was computed for axes ${JSON.stringify(axes)}, but this volume's are ${JSON.stringify(dims.map((d) => d.basis?.vec ?? d.label))}. Recompute it with the volume's nexplan_grid (open_volume).`);
}

export const readPoints = defineTool({
  name: 'read_points',
  title: 'Values at points of a volume',
  description:
    'The processed value at given points of an open volume (the viewer\'s hover readout on a one-bin slice: the equal-weight mean of the distinct measured, unmasked voxels in the symmetry orbit of the voxel containing the point), with the maximum, mean and median over a box of ±box voxels around it. ' +
    'Points are in the volume\'s axis coordinates, or come from a NEXPLAN bragg_positions file (points_path, JSON or CSV, computed for this volume\'s grid), whose hkl and class (present, systematic absence, satellite, …) are kept and summarised: e.g. to see whether forbidden nodes hold intensity. NaN means no data.',
  input: z.strictObject({
    volume_id: volumeId,
    points: z.array(vec3).max(100000).optional().describe('Points in the volume\'s axis coordinates.'),
    points_path: z.string().optional().describe('A NEXPLAN bragg_positions output (.json or .csv) for this volume\'s grid.'),
    box: z.number().int().min(0).max(5).optional().describe('Half-width of the box in voxels (default 1: 3×3×3).'),
    limit: z.number().int().min(1).max(100000).optional().describe('Most points returned (default 200; output_path gets all).'),
    output_path: outputPath.optional().describe('Write every point as .csv.'),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  run(args, { workspace }) {
    if ((args.points === undefined) === (args.points_path === undefined)) throw new ToolError('Give points or points_path.');
    const d = workspace.get(args.volume_id), { dims, shape } = d.info;
    const nodes = args.points ? args.points.map((at) => ({ at })) : readNodes(args.points_path, dims);
    const n = [shape[2], shape[1], shape[0]], strides = [1, n[0], n[0] * n[1]];
    const e0 = dims.map((x) => x.edges[0]), w = dims.map(binWidth);
    const mask = d.mask?.bits ?? null, maps = d.symmetry.maps, volume = d.volume;
    const members = new Set();
    // The symmetrized, masked value of the voxel at bin indices p (NaN without data).
    const orbitMean = (p) => {
      members.clear();
      for (const { M, t } of maps) {
        const q = [0, 1, 2].map((i) => M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + t[i]);
        if (q.every((qi, i) => qi >= 0 && qi < n[i])) members.add(q[0] + q[1] * strides[1] + q[2] * strides[2]);
      }
      let sum = 0, count = 0;
      for (const f of members) {
        const v = volume[f];
        if (v === v && !(mask && mask[f])) { sum += v; count++; }
      }
      return count ? sum / count : NaN;
    };
    const box = args.box ?? 1;
    const rows = nodes.map((node) => {
      const p = node.at.map((x, i) => Math.floor((x - e0[i]) / w[i]));
      const inside = p.every((pi, i) => pi >= 0 && pi < n[i]);
      const row = { ...(node.hkl ? { hkl: node.hkl } : {}), at: node.at, ...(node.class ? { class: node.class } : {}) };
      if (!inside) return { ...row, voxel: null, value: NaN, note: 'outside the volume' };
      const values = [];
      for (let a = -box; a <= box; a++) for (let b = -box; b <= box; b++) for (let c = -box; c <= box; c++) {
        const q = [p[0] + a, p[1] + b, p[2] + c];
        if (q.every((qi, i) => qi >= 0 && qi < n[i])) {
          const v = orbitMean(q);
          if (v === v) values.push(v);
        }
      }
      values.sort((x, y) => x - y);
      return {
        ...row,
        voxel: p,
        value: orbitMean(p),
        box_max: values.length ? values[values.length - 1] : NaN,
        box_mean: values.length ? values.reduce((s, v) => s + v, 0) / values.length : NaN,
        box_median: values.length ? values[Math.floor(values.length / 2)] : NaN,
        box_voxels: values.length,
      };
    });
    const median = (xs) => {
      const s = xs.filter((x) => x === x).sort((x, y) => x - y);
      return s.length ? s[Math.floor(s.length / 2)] : NaN;
    };
    const classes = {};
    for (const r of rows) {
      const k = r.class ?? 'points';
      (classes[k] ??= []).push(r);
    }
    // Medians over the points with data at the node itself.
    const summary = Object.fromEntries(Object.entries(classes).map(([k, list]) => {
      const measured = list.filter((r) => r.value === r.value);
      return [k, {
        points: list.length,
        with_data: measured.length,
        median_value: median(measured.map((r) => r.value)),
        median_box_max: median(measured.map((r) => r.box_max)),
        median_box_median: median(measured.map((r) => r.box_median)),
      }];
    }));
    const limit = args.limit ?? 200;
    const out = {
      volume_id: d.id,
      box: `±${box} voxels`,
      processing: d.processing(),
      summary,
      points: rows.slice(0, limit),
      ...(rows.length > limit ? { truncated: `${rows.length - limit} more points; give output_path for all` } : {}),
    };
    if (args.output_path) {
      const labels = dims.map((x) => x.label);
      const lines = ['h,k,l,' + labels.map((l) => `at_${l}`).join(',') + ',i0,i1,i2,class,value,box_max,box_mean,box_median,box_voxels',
        ...rows.map((r) => [...(r.hkl ?? ['', '', '']), ...r.at.map(num), ...(r.voxel ?? ['', '', '']), r.class ?? '', ...[r.value, r.box_max, r.box_mean, r.box_median].map(num), r.box_voxels ?? 0].join(','))];
      out.written = writeOutput(args.output_path, `${lines.join('\n')}\n`, args.overwrite);
    }
    return out;
  },
});
