// Opening volumes and describing them: the file's axes, cell and uncertainties, and the grid in the
// form NEXPLAN's volume tools take, so Bragg nodes and coverage maps can be computed for it.

import { z } from 'zod';

import { isHKL } from '../../js/nexus.js';
import { binWidth, isUniform, planeName, spaceOf, valueStats, volumeId } from '../shared.js';
import { defineTool } from '../tool.js';

function formatOf(info) {
  if (info.signal.startsWith('/MDHistoWorkspace/')) return 'Mantid MDHistoWorkspace (SaveMD)';
  if (/\/data$/.test(info.signal) && info.dims.some((d) => /^[hkxyz]_axis$/.test(d.name))) {
    return info.dims.every((d) => d.frame === 'direct') ? 'NEBULA3D 3D-ΔPDF' : 'NEBULA3D volume';
  }
  return 'NXdata';
}

/** Mantid MDNorm binning text "min,step,max" for uniform edges. */
const binning = (dim) => {
  const n = (x) => String(Number(x.toPrecision(10)));
  return `${n(dim.edges[0])},${n(binWidth(dim))},${n(dim.edges[dim.edges.length - 1])}`;
};

/** The grid as NEXPLAN's bragg_positions, coverage_map and laue_symmetry take it, or why it cannot be given. */
function nexplanGrid(dims) {
  if (!dims.every((d) => d.basis && isHKL(d))) return { unavailable: 'NEXPLAN grids need axes along hkl directions (r.l.u.).' };
  if (!dims.every(isUniform)) return { unavailable: 'NEXPLAN grids need uniform bins.' };
  return {
    axes: dims.map((d) => d.basis.vec),
    binning: dims.map(binning),
    note: 'Add q_convention: "inelastic" when the volume was reduced with Mantid\'s default Q.convention (most SNS MDNorm data), "crystallography" when it was reduced with the crystallographic one (usual for X-ray data).',
  };
}

export function describe(d, { full = true } = {}) {
  const { info } = d, l = info.lattice;
  const out = {
    volume_id: d.id,
    file: d.path,
    format: formatOf(info),
    space: spaceOf(info.dims),
    shape: info.shape,
    axes: info.dims.map((dim, i) => ({
      index: i,
      label: dim.label,
      long_name: dim.longName,
      units: dim.units,
      min: dim.edges[0],
      max: dim.edges[dim.edges.length - 1],
      bins: dim.edges.length - 1,
      step: binWidth(dim),
      ...(isUniform(dim) ? {} : { uniform: false }),
      ...(dim.basis ? { direction_hkl: dim.basis.vec } : {}),
    })),
    planes: [2, 1, 0].map((fixed) => ({ plane: planeName(info.dims, fixed), fixed_axis: info.dims[fixed].label })),
    measured_fraction: d.stats.fraction,
    processing: d.processing(),
  };
  if (!full) return out;
  return {
    ...out,
    signal: info.signal,
    voxels: info.voxels,
    values: valueStats(d.volume),
    cell: l ? { a: l.a, b: l.b, c: l.c, alpha: l.alpha, beta: l.beta, gamma: l.gamma, source: l.source } : null,
    ...(l?.ub ? { ub: [l.ub.slice(0, 3), l.ub.slice(3, 6), l.ub.slice(6, 9)], ub_convention: 'Mantid: q = UB·h without 2π (NEBULA3D stores 2π·UB)' } : {}),
    uncertainties: info.errors ? `${info.errors.squared ? 'variances' : 'standard deviations'} in ${info.errors.path}` : 'none',
    ...(info.signed ? { signed: 'a ΔPDF: values change sign' } : {}),
    nexplan_grid: nexplanGrid(info.dims),
    storage: { order: 'C order of shape; display axis 0 is the fastest storage axis (Mantid D0), as the viewer and Mantid store it', dtype: info.dtype, chunks: info.chunks },
  };
}

export const openVolume = defineTool({
  name: 'open_volume',
  title: 'Open a 3-D volume',
  description:
    'Opens a 3-D histogram from a NeXus/HDF5 file and keeps it in memory for the other tools: a Mantid MDHistoWorkspace (SaveMD, e.g. MDNorm output), any NXdata group with a 3-D signal, or a NEBULA3D volume or 3D-ΔPDF. ' +
    'Returns its id, axes (labels, ranges, bins, hkl directions), unit cell and UB, value statistics, the slice planes, and the grid in the form NEXPLAN\'s bragg_positions, coverage_map and laue_symmetry take. ' +
    'Masked and non-finite voxels count as unmeasured. Opening another file under an id in use replaces it.',
  input: z.strictObject({
    path: z.string().min(1).describe('The .nxs / .h5 file.'),
    id: z.string().regex(/^[\w.-]{1,64}$/).optional().describe('Id for the volume (default: from the file name).'),
    align_to: z.string().optional().describe('Id of an open volume whose axis order to use, when this file has the same axes in another order (for comparing).'),
  }),
  annotations: { readOnlyHint: true, openWorldHint: false },
  async run(args, { workspace }) {
    const d = await workspace.open(args.path, { id: args.id, alignTo: args.align_to });
    return { ...describe(d), loaded_in_s: d.seconds };
  },
});

export const describeVolume = defineTool({
  name: 'describe_volume',
  title: 'Describe an open volume',
  description: 'An open volume\'s axes, cell, statistics and current processing (symmetry and mask); without volume_id, the volumes open in this session.',
  input: z.strictObject({ volume_id: volumeId.optional() }),
  annotations: { readOnlyHint: true, openWorldHint: false },
  run(args, { workspace }) {
    if (args.volume_id === undefined) return { volumes: workspace.list().map((d) => describe(d, { full: false })) };
    const d = workspace.get(args.volume_id);
    return { ...describe(d), symmetry: d.describeSymmetry() };
  },
});

export const closeVolume = defineTool({
  name: 'close_volume',
  title: 'Close a volume',
  description: 'Frees an open volume\'s memory. Its id can then be reused.',
  input: z.strictObject({ volume_id: volumeId }),
  annotations: { readOnlyHint: true, openWorldHint: false },
  run(args, { workspace }) {
    workspace.close(args.volume_id);
    return { closed: args.volume_id, open: workspace.list().map((d) => d.id) };
  },
});
