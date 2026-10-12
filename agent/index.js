// The NeXus Viewer's agent tools: the registry, a validated call, and the tools as JSON Schema
// definitions for any tool-use API (Anthropic `tools`, OpenAI-style function calling). server.js
// serves the same registry over MCP.

import { z } from 'zod';

import { imagesOf, ToolError } from './tool.js';
import { lineCutTool, powderAverageTool, readPoints, sliceVolume } from './tools/analysis.js';
import { exportNebula3d } from './tools/export.js';
import { buildMask, setSymmetry } from './tools/processing.js';
import { closeVolume, describeVolume, openVolume } from './tools/volumes.js';

export { ToolError, toJson } from './tool.js';
export { Workspace } from './workspace.js';

/** Every tool, in the order a session usually needs them. */
export const TOOLS = [
  openVolume,
  describeVolume,
  setSymmetry,
  buildMask,
  sliceVolume,
  lineCutTool,
  powderAverageTool,
  readPoints,
  exportNebula3d,
  closeVolume,
];

/** How to use the tools, for a client's system prompt (MCP `instructions`). */
export const INSTRUCTIONS = `The NeXus Viewer's analysis of reduced 3-D scattering volumes (neutron or X-ray; reciprocal-space HKL histograms such as Mantid MDNorm output, and NEBULA3D volumes and 3D-ΔPDFs), with the viewer's own code, so results match what the viewer shows.

Workflow: open_volume (a .nxs/.h5 path) returns a volume_id and the file's axes, cell and grid. Processing is set on the volume and applies to every later result, as in the viewer's Processing panel:
- set_symmetry: Laue class preset or hkl operations to pool symmetry-equivalent voxels (NEXPLAN's laue_symmetry gives them for a structure).
- build_mask: remove coverage-edge voxels (edge_radius) and symmetry outliers (outlier_k, after set_symmetry).
Then:
- slice_volume: a 2-D slice (HK, HL, KL, … by the axis held fixed) with statistics and a PNG in true geometry to look at.
- line_cut: a 1-D profile along any line, averaged over a rod, with σ.
- powder_average: I(|Q|) with coverage and σ.
- read_points: values at points, or at the Bragg nodes of a NEXPLAN bragg_positions file (e.g. forbidden nodes, satellites) for this grid; open_volume's nexplan_grid is the grid to give NEXPLAN as it is, adding q_convention as its nexplan_grid_note says.
- export_nebula3d: the processed volume as a NEBULA3D input file for the 3D-ΔPDF.

Conventions: coordinates are the volume's own axes (labels such as H, K, L or [H,H,0]); |Q| in Å⁻¹ with 2π; the UB in Mantid's convention, without 2π (a Mantid file's as stored; for a NEBULA3D file, which stores 2π·UB, its ub_matrix divided by 2π); arrays (masks) in the volume's display order, axis 0 fastest, which is the file's own order unless the volume was opened with align_to (then describe_volume's storage gives the file's order too). Files written need overwrite: true to replace an existing file. Several volumes can be open at once (e.g. temperatures): set the same processing on each before comparing.

Limits: nothing is interpolated or smoothed; values are the file's (no background subtraction or absolute scaling). Masked and non-finite voxels are unmeasured. Symmetry averaging needs uniform bins that the operations map onto themselves. Check anything you publish in the viewer itself.`;

/** Validates the arguments and runs the tool; a failure comes back as a message the caller can act on, not a throw. */
export async function runTool(tool, args, ctx) {
  const parsed = tool.input.safeParse(args ?? {});
  if (!parsed.success) return { ok: false, error: `Invalid arguments for ${tool.name}:\n${z.prettifyError(parsed.error)}` };
  try {
    const value = await tool.run(parsed.data, ctx);
    return { ok: true, value, images: imagesOf(value) };
  } catch (e) {
    if (e instanceof ToolError) return { ok: false, error: e.message, ...(e.details ? { details: e.details } : {}) };
    return { ok: false, error: `${tool.name} failed: ${e instanceof Error ? e.message : String(e)}` };
  }
}

/** Runs a tool by name. */
export function callTool(name, args, ctx) {
  const tool = TOOLS.find((t) => t.name === name);
  if (!tool) return Promise.resolve({ ok: false, error: `No tool '${name}'. Tools: ${TOOLS.map((t) => t.name).join(', ')}.` });
  return runTool(tool, args, ctx);
}

/** The tools as JSON Schema (draft 7) definitions: `{ name, description, input_schema }`, the Anthropic Messages API shape. */
export function toolDefinitions() {
  return TOOLS.map((t) => {
    const { $schema: _schema, ...input_schema } = z.toJSONSchema(t.input, { target: 'draft-7', io: 'input' });
    return { name: t.name, description: t.description, input_schema };
  });
}
