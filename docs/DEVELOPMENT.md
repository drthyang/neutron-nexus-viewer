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
| `js/app.js` | UI: control panel, views, layouts, interaction and canvas drawing |
| `js/worker.js` | Module worker: opens the file with h5wasm, holds the volume, answers slice, mask and isosurface requests |
| `js/nexus.js` | Finds the histogram in the file, reads it, and handles the UB matrix and lattice geometry |
| `js/slab.js` | Slab selection and symmetry-pooled slab averaging (pure functions) |
| `js/symmetry.js` | Parsing operations, group closure, metric check and integer index maps |
| `js/mask.js` | Coverage-edge erosion and symmetry-outlier masks |
| `js/iso.js` | Origin-aligned coarse binning, orbit means and surface nets |
| `js/view3d.js` | three.js scene, loaded when the 3-D view opens |
| `js/colormaps.js` | Colormap lookup tables (generated) |

## Fixtures, example data and generated files

These scripts need numpy and h5py, and matplotlib for the colormaps:

- `tools/make_fixtures.py` writes `tests/fixtures/` and their expected slices, computed by a verbatim copy of the reference `average_slab()`.
- `tools/make_example.py` writes `examples/demo_hexagonal.nxs`, a 101³ synthetic hexagonal dataset (1.3 MB) with coverage wedges and bright edge voxels.
- `tools/make_colormaps.py` writes `js/colormaps.js` from matplotlib.

`docs/screenshot.png` is a capture of the example data with 6/mmm averaging and a mask (edge 1, 5σ).

## Dependencies

- [h5wasm](https://github.com/usnistgov/h5wasm) (NIST) reads HDF5 in the worker. It is loaded from jsDelivr at the version pinned in `js/worker.js`, which must match `package.json`, where the tests use it.
- [three.js](https://threejs.org) renders the 3-D view. It is loaded from jsDelivr through the import map in `index.html`.

No other runtime dependencies. The file never leaves the browser; only these two libraries are fetched from the CDN.

## Deployment

GitHub Pages serves the root of `main` (*Deploy from a branch*, `main`, `/ (root)`) at [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/). Every push to `main` redeploys it within about a minute. `.nojekyll` makes Pages serve the files as they are.
