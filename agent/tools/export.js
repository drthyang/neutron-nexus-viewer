// The hand-off to NEBULA3D: the masked, symmetrized volume on a grid symmetric about Q = 0, written as
// the HDF5 file NEBULA3D's 3D-ΔPDF pipeline reads, exactly as the viewer's Export section writes it.

import { existsSync, mkdirSync } from 'node:fs';
import { dirname } from 'node:path';

import h5wasm from 'h5wasm/node';
import { z } from 'zod';

import { exportPlan, sigmaNote, symmetrizeForExport, writeNebulaFile } from '../../js/export.js';
import { formatOp } from '../../js/symmetry.js';
import { resolvePath } from '../files.js';
import { outputPath, overwrite, volumeId, WRITTEN_BY } from '../shared.js';
import { defineTool, ToolError } from '../tool.js';

export const exportNebula3d = defineTool({
  name: 'export_nebula3d',
  title: 'Export a volume for NEBULA3D',
  description:
    'Writes an open volume, masked and symmetrized as set (set_symmetry, build_mask), as a NEBULA3D input file: a Mantid MDHistoWorkspace (SaveMD layout, which Mantid and other NeXus readers open too) as NEBULA3D writes its own, in (H, K, L) C order (D2 = H, D1 = K, D0 = L), on a grid padded to be symmetric about 0 (NEBULA3D takes the centre of the middle bin as Q = 0), with the UB in the oriented lattice and NEBULA3D\'s nebula3d group holding the exact axes and UB (with 2π). ' +
    'Every voxel of the padded grid gets the mean of the measured, unmasked voxels among its symmetry equivalents, so symmetry also fills the padding; voxels with none are written as 0 and masked (Mantid\'s mask: 1 = masked) for NEBULA3D to backfill. errors_squared is that mean\'s variance, (Σσ²)/n² from the file\'s uncertainties; a file without them gets |data|, the √|data| NEBULA3D takes as σ, marked sigma = sqrt_intensity. The source file, symmetry operations and mask are recorded as attributes of the nebula3d group. Needs plain H, K, L axes with uniform bins, a bin centred at 0 on each, and a unit cell; a grid with a bin edge at 0 is refused with the binning to reduce with instead (rebinning here would smooth the data).',
  input: z.strictObject({
    volume_id: volumeId,
    output_path: outputPath.describe('The .nxs / .h5 file to write.'),
    overwrite,
  }),
  annotations: { readOnlyHint: false, destructiveHint: false, openWorldHint: false },
  async run(args, { workspace }) {
    const d = workspace.get(args.volume_id), { dims, lattice, shape } = d.info;
    let plan;
    try {
      plan = exportPlan(dims, lattice);
    } catch (err) {
      throw new ToolError(err.message);
    }
    const file = resolvePath(args.output_path);
    if (existsSync(file) && !args.overwrite) throw new ToolError(`${file} exists. Pass overwrite: true to replace it, or another output_path.`);
    const maps = d.mapsOn(plan.paddedDims);
    const t0 = performance.now();
    const variance = d.readVariance();
    const result = symmetrizeForExport(d.volume, shape, plan, maps, { mask: d.mask?.bits ?? null, variance });
    const sym = d.symmetry.ops.length > 1;
    const uncertainties = sigmaNote(result, d.info.errors?.path.split('/').pop(), d.varianceNote);
    const attrs = {
      source_file: d.name,
      symmetry: sym ? d.symmetry.name : 'none',
      symmetry_ops: d.symmetry.ops.map(formatOp).join('; '),
      mask: d.processing().mask,
      ub_source: lattice.source,
      padded: plan.padded ? 'yes' : 'no',
      uncertainties,
      created_by: WRITTEN_BY,
      created: new Date().toISOString(),
    };
    await h5wasm.ready;
    mkdirSync(dirname(file), { recursive: true });
    writeNebulaFile(h5wasm, file, plan, result, attrs);
    const { stats } = result;
    return {
      volume_id: d.id,
      written: { path: file },
      shape_hkl: stats.shape,
      padded: plan.padded,
      valid_voxels: stats.valid,
      valid_fraction: stats.valid / stats.total,
      measured_voxels_used: stats.measured,
      sigma: !!result.sigma,
      attributes: attrs,
      seconds: (performance.now() - t0) / 1000,
      next: 'Open the file in NEBULA3D (its web console\'s file picker, or run_pipeline on it).',
    };
  },
});
