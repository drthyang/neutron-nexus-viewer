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
| `index.html` | Markup and styles (light variant of the nebula3d design tokens) |
| `js/app.js` | UI: control panel, views, layouts, interaction, canvas drawing and the second (comparison) dataset |
| `js/worker.js` | Module worker (one per open dataset): opens the file with h5wasm, holds the volume, answers slice, mask and isosurface requests |
| `js/nexus.js` | Finds the histogram in the file, reads it, and handles the UB matrix and lattice geometry |
| `js/slab.js` | Slab selection and symmetry-pooled slab averaging (pure functions) |
| `js/symmetry.js` | Parsing operations, group closure, metric check and integer index maps |
| `js/mask.js` | Coverage-edge erosion and symmetry-outlier masks |
| `js/iso.js` | Origin-aligned coarse binning, orbit means and surface nets |
| `js/export.js` | NEBULA3D export: grid plan and padding, full-volume symmetrization, HDF5 writer |
| `js/view3d.js` | three.js scene, loaded when the 3-D view opens |
| `js/colormaps.js` | Colormap lookup tables (generated) |

## Fixtures, example data and generated files

These scripts need numpy and h5py, and matplotlib for the colormaps:

- `tools/make_fixtures.py` writes `tests/fixtures/` and their expected slices, computed by a verbatim copy of the reference `average_slab()`.
- `tools/make_example.py` writes the example: `examples/demo_300K.nxs` and `examples/demo_10K.nxs`, one 101³ synthetic hexagonal crystal (about 1.5 MB each) above and below an ordering transition. Both have Bragg peaks (with thermal diffuse halos at 300 K), M-point short-range-order rods (300 K) or superlattice peaks (10 K), coverage wedges and bright edge voxels.
- `tools/make_colormaps.py` writes `js/colormaps.js` from matplotlib.

`docs/screenshot.png` is a capture of the example (`?demo`: 6/mmm averaging and a 1-voxel edge mask) at 1440×900 and 1.5× pixel ratio.

## Dependencies

- [h5wasm](https://github.com/usnistgov/h5wasm) (NIST) reads HDF5 in the worker. It is loaded from jsDelivr at the version pinned in `js/worker.js`, which must match `package.json`, where the tests use it.
- [three.js](https://threejs.org) renders the 3-D view. It is loaded from jsDelivr through the import map in `index.html`.

No other runtime dependencies. The file never leaves the browser; only these two libraries are fetched from the CDN.

## Deployment

GitHub Pages serves the root of `main` (*Deploy from a branch*, `main`, `/ (root)`) at [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/). Every push to `main` redeploys it within about a minute. `.nojekyll` makes Pages serve the files as they are.

Pages lets browsers cache files for 10 minutes. `index.html` therefore loads `js/app.js` with its own `Last-Modified` date as a query string, so a freshly deployed page never runs with an older cached `app.js`. The other modules are loaded by plain relative URLs.
