# Agent tools and the MCP server

[← README](../README.md) · [Features](FEATURES.md) · [User guide](USER_GUIDE.md) · [Method](METHOD.md) · [Development](DEVELOPMENT.md)

The viewer's analysis is also available as tools for AI agents. The tools are in `agent/`, and
`agent/main.js` serves them as a [Model Context Protocol](https://modelcontextprotocol.io) server on
stdio, so Claude Code, Claude Desktop or any other MCP client can use them. An agent can open a
reduced volume, set symmetry averaging and a mask, look at slices, take line cuts and I(Q), read the
intensity at Bragg nodes from NEXPLAN, and export the processed volume for NEBULA3D, all by calling
tools.

The tools call the viewer's own modules (`js/nexus.js`, `slab.js`, `cut.js`, `symmetry.js`,
`mask.js`, `powder.js`, `export.js`), so for the same choices they give the viewer's numbers. The
page itself does not change and still needs nothing installed.

## Setup

Install the dependencies once: `npm install` (Node 20.10 or later).

**Claude Code, in this repository.** `.mcp.json` registers the server as `nexus-viewer`; Claude Code
asks once whether to enable it.

**Claude Code, in another project, or Claude Desktop.** Register `agent/main.js` with its absolute
path:

```bash
claude mcp add --scope user nexus-viewer -- node /path/to/neutron-nexus-viewer/agent/main.js
```

```json
{ "mcpServers": { "nexus-viewer": { "command": "node", "args": ["/path/to/neutron-nexus-viewer/agent/main.js"] } } }
```

There is no build step: the server runs the files as they are, so it always matches the code.

**With the other apps.** The [scattering-agent-tools](https://github.com/drthyang/scattering-agent-tools)
plugin starts this server from this checkout, next to NEXPLAN's and MATERIA's, and adds the skills
that combine them.

## Tools

| Tool | What it does | In the viewer |
| --- | --- | --- |
| `open_volume` | Opens a Mantid MDHistoWorkspace, NXdata or NEBULA3D file; returns its id, axes, cell and UB, value statistics, slice planes, and the grid in NEXPLAN's form (`nexplan_grid`). | Open a file |
| `describe_volume` | One volume's description and processing, or the volumes open. | File details |
| `set_symmetry` | Laue class preset or hkl operations to pool equivalent voxels; fails when they do not fit the grid, warns when they change the metric. | Processing → Symmetry averaging |
| `build_mask` | Coverage-edge erosion and symmetry-outlier rejection; writes the mask as `.npy`. | Processing → Mask |
| `slice_volume` | A slice by the axis held fixed: statistics, Auto color range, a PNG in true geometry; values as `.npy`, `.csv` or `.png`. | Slice views |
| `line_cut` | A 1-D profile along any line, averaged over a rod, with σ; `.csv` or the viewer's text format. | Cut mode |
| `powder_average` | I(\|Q\|) with σ and shell coverage; `.csv` or the viewer's text format. | I(Q) |
| `read_points` | Values at points, with box statistics; takes NEXPLAN `bragg_positions` files and summarises by class (present, systematic absence, satellite, …). | Hover readout |
| `export_nebula3d` | The masked, symmetrized volume, with σ, as a NEBULA3D input file. | Export → Open in NEBULA3D |
| `close_volume` | Frees a volume's memory. | — |

The server also gives `docs/METHOD.md` as the resource `nexus-viewer://docs/method`, and its
`instructions` (`INSTRUCTIONS` in `agent/index.js`) tell the client how the tools fit together.

## How the tools work

- **A session.** `open_volume` keeps the volume in memory and returns its id (from the file name, or
  `id`); the other tools take `volume_id`. Several volumes can be open, e.g. two temperatures;
  `align_to` shows a file's axes in another volume's order.
- **Processing.** `set_symmetry` and `build_mask` stay set on the volume and apply to every later
  slice, cut, I(Q), point readout and export, as the viewer's Processing panel does. Each result
  says which processing it used. Comparing two volumes needs the same processing set on both.
- **Coordinates.** Points and slab positions are in the volume's own axes (their labels: H, K, L, or
  e.g. [H,H,0]). Line-cut widths are in Å⁻¹ with 2π when the file has a cell, as in the viewer.
  Typed positions are placed on voxels: `line_cut` slides its ends along the line to voxel centres
  of its axis (and says so), the default slab is the bin at 0, and `read_points` flags a point on a
  bin edge (an integer node on a grid with edges at the integers) and averages the voxels that meet
  there, with its box centred on the point.
  `angles` chooses nominal (default, as the viewer's switch) or measured cell angles for drawing and
  distances.
- **Images.** `slice_volume` returns a PNG (unless `image: false`): the slice drawn as the viewer draws
  it, with oblique axes at their real angle, x to the right and y up, no smoothing, grey where there
  is no data, and the viewer's color scales (Auto: vmin 0, vmax at the 97th percentile, asinh
  softening at the median).
- **Files.** Inputs are read from the local disk. `output_path` writes a file; an existing file is an
  error unless `overwrite: true`. Relative paths start from the directory the client started the
  server in; `~` is the home directory. Text files start with `#` lines saying what was computed and
  with which processing.
- **Results.** JSON, numbers to 7 significant figures, NaN as `null`. Long lists are cut to `limit`
  entries (the file gets all). A failure is an MCP tool error whose text says what to change, with
  the details (open ids, presets) as JSON. Unknown arguments are refused.

## With NEXPLAN and NEBULA3D

- `open_volume`'s `nexplan_grid` is the `grid` NEXPLAN's `laue_symmetry`, `bragg_positions` and
  `coverage_map` take; add `q_convention` (`inelastic` for Mantid's default, `crystallography`
  otherwise).
- NEXPLAN's `laue_symmetry` gives `nexus_viewer_preset` (for `set_symmetry`'s `laue_class`) or
  `symmetry_ops` (its `operations`).
- `read_points` with `points_path` reads NEXPLAN's `bragg_positions` output (JSON or CSV) for the same
  grid and summarises the values by node class: intensity at systematic absences points to λ/2,
  multiple scattering, a lower symmetry or a superlattice.
- `export_nebula3d` writes the file NEBULA3D's pipeline reads, with σ propagated from the file's
  uncertainties in `/entry/sigma` and the symmetry operations in `/entry@symmetry_ops`.

## Using the tools from code

`agent/index.js` exports the registry without a transport:

- `TOOLS`: each tool's `name`, `description`, zod `input` schema and `run`.
- `toolDefinitions()`: `{ name, description, input_schema }` with JSON Schema (draft 7), the shape of
  the Claude Messages API `tools` parameter.
- `callTool(name, input, { workspace })`: validates the input and runs the tool; it returns
  `{ ok: true, value, images }` or `{ ok: false, error, details }` and never throws for a tool failure.

To add a tool: write it with `defineTool` in `agent/tools/`, add it to `TOOLS`, document it in the
table above, and test it in `tests/agent.test.js` (the tests check that every tool is documented
here).

## Limits

The tools have the viewer's limits: values are the file's, without background subtraction or
absolute scaling; nothing is interpolated; symmetry averaging needs uniform bins that the operations
map onto themselves; I(Q) needs a cell (or Q axes); the export needs plain H, K, L axes with a bin centred at 0. The Rigaku
reduction, reindexing and the 3-D view are in the page only for now. An agent can still misread a
result, so check what you publish in the viewer.
