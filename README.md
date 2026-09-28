# NeXus Viewer

[![Live demo](https://img.shields.io/badge/demo-live-2563eb)](https://drthyang.github.io/neutron-nexus-viewer/)
[![Tests](https://github.com/drthyang/neutron-nexus-viewer/actions/workflows/test.yml/badge.svg)](https://github.com/drthyang/neutron-nexus-viewer/actions/workflows/test.yml)
[![Runs in your browser](https://img.shields.io/badge/runs%20in-your%20browser-6b6b6b)](https://drthyang.github.io/neutron-nexus-viewer/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)

**Reciprocal-space slices of 3-D neutron histograms, with symmetry averaging, artifact masking and a 3-D view — all in your browser.**

**▶ Try it: [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/)** — open a Mantid `.nxs` file or the [example data](https://drthyang.github.io/neutron-nexus-viewer/?url=examples/demo_hexagonal.nxs). Nothing to install, and your file never leaves your machine.

<p align="center"><img src="docs/screenshot.png" alt="NeXus Viewer with a hexagonal example dataset: the control panel shows the processing pipeline with 6/mmm symmetry averaging and a mask that removed 19% of voxels; the workspace shows the HK, HL and KL slices and a 3-D isosurface with the three slice planes." width="100%"/></p>

## Features

- **Orthogonal slices** — HK, HL and KL cuts with live slab center and thickness, drawn in the true lattice geometry from the UB matrix. Click to move the other two slices through a point, drag a box to zoom, or drag to pan.
- **Symmetry averaging** — any Laue class or your own operations, closed into a group and applied exactly on the bin grid, with no interpolation. A warning appears when the operations do not fit the cell.
- **Artifact masking** — removes detector-edge voxels and symmetry outliers before averaging. Preview what is removed, and export the mask as a NumPy array.
- **3-D view** — a transparent isosurface of the processed volume, with the current slices as planes.
- **Private and fast** — the file is read locally with [h5wasm](https://github.com/usnistgov/h5wasm). A 401³ volume loads in about 2 s, and slices update in tens of milliseconds.

Reads Mantid `MDHistoWorkspace` files written by `SaveMD`, and any NeXus/HDF5 file whose `NXdata` group holds a 3-D histogram. The [user guide](docs/USER_GUIDE.md#supported-files) lists the details.

## Run locally

The app is a static page with no build step:

```bash
git clone https://github.com/drthyang/neutron-nexus-viewer.git
cd neutron-nexus-viewer
python3 -m http.server 8000   # open http://localhost:8000
```

Tests (`npm install && npm test`), code layout and deployment are covered in [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md).

## Documentation

| Document | What it covers |
| --- | --- |
| [User guide](docs/USER_GUIDE.md) | The interface, every control, layouts, export, remote files and supported formats |
| [Method](docs/METHOD.md) | How slices, symmetry averaging, masking and the 3-D view are computed, and how they are validated |
| [Development](docs/DEVELOPMENT.md) | Code layout, tests and fixtures, example data, dependencies and deployment |

## License

Released under the [GNU Affero General Public License v3.0](LICENSE) © 2026 Tsung-Han Yang.

The AGPL is a strong copyleft license: you may use, study, modify, and redistribute this
software, but derivative works must also be released under the AGPLv3. Notably, if you run a
modified version as a **network service**, you must offer its complete source code to the users
of that service (AGPL §13).

*This project is personal work, developed and maintained in my personal capacity.*
