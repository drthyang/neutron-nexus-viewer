# NeXus Viewer

[![Live demo](https://img.shields.io/badge/demo-live-2563eb)](https://drthyang.github.io/neutron-nexus-viewer/)
[![Tests](https://github.com/drthyang/neutron-nexus-viewer/actions/workflows/test.yml/badge.svg)](https://github.com/drthyang/neutron-nexus-viewer/actions/workflows/test.yml)
[![License: AGPL v3](https://img.shields.io/badge/license-AGPL--3.0-3c8c3c)](LICENSE)

A browser app for exploring 3-D neutron and X-ray scattering volumes: reciprocal-space slices,
line cuts and I(Q) of reduced data, with nothing to install and files that never leave your machine.

<p align="center">
  <img src="docs/screenshot.png" alt="NeXus Viewer comparing the example, a synthetic hexagonal crystal, at 300 K and 10 K: the HK, HL and KL slices are each split along their diagonal, with 300 K in the lower-left half and 10 K in the upper-right half. A line cut drawn on the HL slice along L through the M point (0.5, 0, L) is plotted in the fourth view: at 300 K a broad, continuous diffuse rod, at 10 K sharp superlattice peaks at even L. The control panel shows the processing pipeline, one line per stage, with 6/mmm symmetry averaging, a mask and the Reindex card, and the header row the click modes, the A / Split / B switch and the color scale. © 2026 Tsung-Han Yang in the bottom-left corner." width="100%" />
</p>

## What it does

- **Slices and line cuts** of Mantid `MDHistoWorkspace`, `NXdata` and NEBULA3D files, drawn in the
  true lattice geometry, with nothing smoothed or interpolated.
- **Symmetry averaging and artifact masks**, applied exactly on the bin grid, with a check that the
  operations fit the cell.
- **Up to four datasets compared** (temperatures, or X-ray next to neutron), with I(Q) and a 3-D view.
- **UB refinement** against the volume's Bragg peaks, with the cell constrained to a crystal system,
  **reindexing with a new UB** (a resampling, which smooths), and **reduction of raw Rigaku XRD
  frames** to an HKL volume.
- **Hand-off to [NEBULA3D](https://github.com/drthyang/nebula3d)** for a 3D-ΔPDF.

## Quick start

Open **[drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/)**
and drop a `.nxs` file on the page, or open [the example](https://drthyang.github.io/neutron-nexus-viewer/?demo),
a synthetic hexagonal crystal at 300 K and 10 K.
To run it locally, serve the folder (`python3 -m http.server 8000`). A first session, step by step,
is in [Features](docs/FEATURES.md#a-first-session).

## Agent tools

`npm install && npm run mcp` serves the analysis to AI agents as 10 MCP tools (Claude Code, Claude
Desktop, any MCP client), computed by the viewer's own code. The tool list and setup are in
[docs/AGENT.md](docs/AGENT.md); [scattering-agent-tools](https://github.com/drthyang/scattering-agent-tools)
installs them together with the other apps' tools.

## Status

No tagged releases yet: the live app is deployed from `main` on every push, and CI runs the tests,
which check the slice engine against the reference Python implementation, on every push.

## Documentation

- [Features](docs/FEATURES.md): the goals, everything the viewer does, and a first session
- [User guide](docs/USER_GUIDE.md): every control, supported files, comparison, reindexing, Rigaku reduction, export
- [Method](docs/METHOD.md): how each result is computed, and how it is validated
- [Agent tools](docs/AGENT.md): the MCP server, its tools and their conventions
- [Development](docs/DEVELOPMENT.md): code layout, tests and fixtures, dependencies, deployment

## License

[GNU Affero General Public License v3.0](LICENSE) © 2026 Tsung-Han Yang. If you run a modified
version on a network server, you must make its complete source available to its users (AGPL §13).

*This project is personal work, developed and maintained in my personal capacity.*
