// The viewer's Processing panel: symmetry averaging and the artifact mask. Both stay set on the volume
// and apply to every later slice, cut, I(Q), point readout and export, as in the viewer.

import { z } from 'zod';

import { edgeMask, maskStats, outlierMask } from '../../js/mask.js';
import { PRESETS } from '../../js/symmetry.js';
import { npyBytes, writeOutput } from '../files.js';
import { outputPath, overwrite, volumeId } from '../shared.js';
import { defineTool, ToolError } from '../tool.js';
import { maskSummary } from '../workspace.js';

const PRESET_NAMES = PRESETS.map(([name]) => name);

export const setSymmetry = defineTool({
  name: 'set_symmetry',
  title: 'Set the symmetry averaging',
  description:
    'Sets the Laue group whose symmetry-equivalent voxels are pooled (equal-weight mean of the distinct measured voxels of each orbit, mapped bin onto bin, no interpolation) in every later slice, line cut, I(Q), point readout, outlier mask and NEBULA3D export of this volume. ' +
    `Give a Laue class preset (${PRESET_NAMES.join(', ')}; trigonal and hexagonal on hexagonal axes) or operations on hkl such as "h+k,-h,l; -k,-h,-l; -h,-k,-l" (generators are closed into a group; x,y,z triplets are converted). ` +
    'NEXPLAN\'s laue_symmetry gives the preset or operations for a structure. Fails when an operation does not map the bin grid onto itself, and warns when the operations change the cell metric by more than 2 % (a wrong setting). "1" or "none" turns averaging off.',
  input: z.strictObject({
    volume_id: volumeId,
    laue_class: z.string().optional().describe(`A preset: ${PRESET_NAMES.join(', ')}.`),
    operations: z.string().optional().describe('Operations on hkl, separated by ";" (e.g. "k,-h+k,l; -h,-k,-l"); "none" for no averaging.'),
  }),
  annotations: { readOnlyHint: true, idempotentHint: true, openWorldHint: false },
  run(args, { workspace }) {
    if ((args.laue_class === undefined) === (args.operations === undefined)) throw new ToolError('Give laue_class or operations.', { presets: PRESET_NAMES });
    const d = workspace.get(args.volume_id);
    const laue = args.laue_class?.trim().toLowerCase() === 'none' ? '1' : args.laue_class;
    const out = d.setSymmetry({ preset: laue, operations: args.operations });
    return {
      volume_id: d.id,
      ...out,
      ...(d.mask?.k ? { mask_note: `The mask's outlier cut was built with ${d.mask.group}; call build_mask again to rebuild it with ${out.group}.` } : {}),
    };
  },
});

export const buildMask = defineTool({
  name: 'build_mask',
  title: 'Mask detector-edge artifacts',
  description:
    'Builds the viewer\'s artifact mask on the unsymmetrized volume, before averaging: removes measured voxels within edge_radius voxels of unmeasured ones (the coverage edges, where detector edges and weak normalization leave bright rims), and voxels more than outlier_k σ above the median of their symmetry equivalents (needs a group of at least 3 operations from set_symmetry). ' +
      'σ is the robust spread 1.4826·MAD, but at least the orbit\'s typical counting error from the file\'s uncertainties (or √counts for integer counts without them), so counting noise in low-count data, where MAD is often 0, is not cut. ' +
    'The volume\'s values are never changed; the mask applies to every later slice, cut, I(Q), point readout and export. Both 0 clears it. Writes the mask as .npy (uint8 in storage order, 1 = coverage edge, 2 = outlier) with output_path.',
  input: z.strictObject({
    volume_id: volumeId,
    edge_radius: z.number().int().min(0).max(20).optional().describe('Erosion radius in voxels (default 0: none; 1 or 2 is typical).'),
    outlier_k: z.number().min(0).optional().describe('Outlier cut in robust σ (default 0: none; 5 is a cautious start).'),
    output_path: outputPath.optional(),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  run(args, { workspace }) {
    const d = workspace.get(args.volume_id);
    const radius = args.edge_radius ?? 0, k = args.outlier_k ?? 0;
    if (k > 0 && d.symmetry.maps.length < 3) {
      throw new ToolError('The outlier cut compares symmetry equivalents: set a group with at least 3 operations first (set_symmetry, e.g. 6/mmm or mmm).');
    }
    const t0 = performance.now();
    if (!(radius > 0 || k > 0)) {
      d.mask = null;
      return { volume_id: d.id, mask: 'none', processing: d.processing() };
    }
    const bits = edgeMask(d.volume, d.info.shape, radius);
    const variance = k > 0 ? d.readVariance() : null;
    if (k > 0) outlierMask(d.volume, d.info.shape, d.symmetry.maps, k, bits, { variance });
    const stats = maskStats(d.volume, bits);
    d.mask = { bits, radius, k, group: d.symmetry.name, stats };
    const out = {
      volume_id: d.id,
      mask: maskSummary(d.mask),
      measured_voxels: stats.measured,
      coverage_edge_voxels: stats.edge,
      outlier_voxels: stats.outlier,
      removed_fraction: stats.measured ? (stats.edge + stats.outlier) / stats.measured : 0,
      ...(k > 0 ? {
        outlier_sigma: variance
          ? 'max(1.4826·MAD, √median σ² of the orbit), σ² from the file; a voxel\'s own σ where both are 0'
          : `1.4826·MAD; for orbits of integer counts at least √median, or the voxel's own √count where both are 0; orbits of other values whose MAD is 0 are not cut (no uncertainties${d.varianceNote ? `: ${d.varianceNote}` : ' in the file'})`,
      } : {}),
      seconds: (performance.now() - t0) / 1000,
      processing: d.processing(),
    };
    if (args.output_path) {
      out.written = writeOutput(args.output_path, npyBytes(bits, d.info.shape), args.overwrite);
      out.array = { dtype: 'uint8', shape: d.info.shape, order: 'storage (Mantid) order, as the volume; 1 = coverage edge, 2 = symmetry outlier' };
    }
    return out;
  },
});
