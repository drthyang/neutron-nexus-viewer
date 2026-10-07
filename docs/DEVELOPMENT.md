# Development

[← README](../README.md) · [User guide](USER_GUIDE.md) · [Method](METHOD.md)

## Run and test

There is no build step. Serve the folder and open it:

```bash
python3 -m http.server 8000   # open http://localhost:8000
```

The tests run in Node, against small HDF5 fixtures:

```bash
npm install && npm test
```

CI runs the same tests on every push (`.github/workflows/test.yml`).

## Code layout

| File | Role |
| --- | --- |
| `index.html` | Markup and styles (light variant of the nebula3d design tokens). Sizes are in rem, 10 px at 1×; a script in its head scales the root size on large screens at 100%, and the canvases draw at the same scale |
| `js/app.js` | UI: control panel, views, layouts, interaction, canvas drawing (the slices, and the 1-D plots of I(Q) and the line cut) and the second (comparison) dataset |
| `js/worker.js` | Module worker (one per open dataset): opens the file with h5wasm, holds the volume, answers slice, line cut, mask, isosurface, I(Q) and export requests |
| `js/nexus.js` | Finds the histogram (and its uncertainties) in the file, reads it, and handles the UB matrix and lattice geometry |
| `js/slab.js` | Slab selection and symmetry-pooled slab averaging (pure functions) |
| `js/cut.js` | Line cuts: the points and band of a cut in a slice, and its symmetry-pooled profile with σ (pure functions) |
| `js/symmetry.js` | Parsing operations, group closure, metric check and integer index maps |
| `js/mask.js` | Coverage-edge erosion and symmetry-outlier masks |
| `js/iso.js` | Origin-aligned coarse binning, orbit means and surface nets |
| `js/powder.js` | I(Q): the \|Q\| metric, shells, and the normalized spherical average with split voxels and propagated uncertainties |
| `js/export.js` | NEBULA3D export: grid plan and padding, full-volume symmetrization, HDF5 writer |
| `js/rigaku-format.js` | Rigaku Oxford Diffraction frames: header, TY6 decoding with the header-statistics check, and the CrysAlisPro text files (.par, crystal.ini, datacoll.ini) |
| `js/rigaku-geometry.js` | Kappa-goniometer and detector geometry, Ewald prediction, cells, and the Levenberg-Marquardt fit |
| `js/rigaku-reduce.js` | The reduction: detector mask, 3-D peak search, geometry refinement, gridding, and the Mantid `SaveMD`-layout writer |
| `js/rigaku-worker.js` | Module worker that reduces an experiment folder and returns the .nxs file |
| `js/rigaku-ui.js` | The *Reduce Rigaku XRD* dialog |
| `js/view3d.js` | three.js scene, loaded when the 3-D view opens |
| `js/colormaps.js` | Colormap lookup tables (generated) |

## Fixtures, example data and generated files

These scripts need numpy and h5py, and matplotlib for the colormaps:

- `tools/make_fixtures.py` writes `tests/fixtures/` (a Mantid `MDHistoWorkspace`, a plain `NXdata`, and a NEBULA3D volume and 3D-ΔPDF) and their expected slices, computed by a verbatim copy of the reference `average_slab()`.
- `tools/make_example.py` writes the example: `examples/demo_300K.nxs` and `examples/demo_10K.nxs`, one 101³ synthetic hexagonal crystal (about 1.5 MB each) above and below an ordering transition. Both have Bragg peaks (with thermal diffuse halos at 300 K), M-point short-range-order rods (300 K) or superlattice peaks (10 K), coverage wedges and bright edge voxels.
- `tools/make_colormaps.py` writes `js/colormaps.js` from matplotlib.

`tests/rigaku-synth.js` writes synthetic Rigaku frames (a TY6 encoder and a frame writer) and simulates a rotation experiment for `tests/rigaku.test.js`; no real frames are in the repository. `tests/rigaku-local.test.js` decodes every frame of a real experiment and checks it against the header statistics when `RIGAKU_DIR` points at one (`RIGAKU_DIR=/path/to/experiment npm test`); without it, it is skipped.

`docs/screenshot.png` is a capture of the example (`?demo`: 6/mmm averaging and a 1-voxel edge mask) at 1440×900 and 1.5× pixel ratio.

## Dependencies

- [h5wasm](https://github.com/usnistgov/h5wasm) (NIST) reads HDF5 in the worker. It is loaded from jsDelivr at the version pinned in `js/worker.js`, which must match `package.json`, where the tests use it.
- [three.js](https://threejs.org) renders the 3-D view. It is loaded from jsDelivr through the import map in `index.html`.

No other runtime dependencies. The file never leaves the browser; only these two libraries are fetched from the CDN.

## Deployment

GitHub Pages serves the root of `main` (*Deploy from a branch*, `main`, `/ (root)`) at [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/). Every push to `main` redeploys it within about a minute. `.nojekyll` makes Pages serve the files as they are.

**Open in NEBULA3D** posts the exported file to a NEBULA3D tab; the protocol is described in NEBULA3D's `docs/web.md` (*Import from the NeXus Viewer*). When NEBULA3D has this page's origin (both on `drthyang.github.io`), the tab opens with `noopener` and the messages go over a `BroadcastChannel`, so the two tabs run in separate browser processes. Same-site tabs that hold a window reference to each other share a renderer process, and reloading or closing the viewer could then also end a NEBULA3D run. On another origin the viewer keeps the window reference and uses `postMessage`. To test against a local NEBULA3D (`cd web && npm run dev:pyodide`), open the viewer with `?nebula3d=http://localhost:5173/` (the `postMessage` path); NEBULA3D accepts localhost senders only in its dev server. To test the `BroadcastChannel` path, serve NEBULA3D's Pages build (`npm run build:pages`, `web/dist` as `/nebula3d/`) and this repository from one local server, and open `?nebula3d=/nebula3d/`.

Pages lets browsers cache files for 10 minutes. `index.html` therefore loads `js/app.js` with its own `Last-Modified` date as a query string, so a freshly deployed page never runs with an older cached `app.js`. The other modules are loaded by plain relative URLs.
