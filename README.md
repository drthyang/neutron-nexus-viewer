# NeXus slice viewer

A static web page for browsing reciprocal-space slices of 3-D neutron histograms, such as Mantid `MDHistoWorkspace` files written by `SaveMD`. Open or drop a `.nxs` file and the page shows three orthogonal slices (HK, HL, KL) with adjustable slab center and thickness, a 3-D view with a transparent isosurface and the current slices, and symmetry averaging with any Laue class or your own operations.

The file is read in your browser by [h5wasm](https://github.com/usnistgov/h5wasm) (HDF5 compiled to WebAssembly) and is never uploaded. The page is static, so it can be hosted on GitHub Pages.

## Using it

- **Open a file** with the button or drag it onto the page. The histogram is loaded into memory as float32: a 401³ volume needs about 260 MB and takes 1–2 s to read.
- **Slices**: each panel has a center and full-thickness control for its integrated axis. Bins whose centers lie within thickness/2 of the center are averaged, exactly as in the original Python viewer. Slices update live while you drag.
- **Click a plot** to move the other two slices through that point. Dashed guides show where they cut.
- **Color**: colormap, asinh, linear or log scale, vmin/vmax, and asinh softening. *Auto range* sets vmax to the 97th percentile and softening to the median of the positive values in the current slices.
- **Symmetry**: pick a Laue class or type generators (see below). Slices and the 3-D view pool every voxel with its symmetry equivalents. The default is no symmetry, so the data are shown as measured.
- **3-D view**: a transparent isosurface of the binned volume with the three current slices as planes, drawn in the lattice geometry and clipped to the view range. Set the level (typed or on a log slider), opacity, and grid (about 64, 100 or 150 blocks per axis); drag to rotate, scroll to zoom, right-drag to pan.
- **Cell angles**: axes are drawn with the reciprocal metric from the UB matrix in the file (or its stored unit cell when there is no UB). *Nominal* snaps direct-cell angles within 1° of 60°, 90° or 120° (for example, 90/90/120 for a hexagonal cell). *Measured* uses them as derived.
- **Save PNG** exports a panel at 3× resolution with its colorbar, without guides.
- **Remote files**: `?url=https://…/file.nxs` downloads and opens a file. The host must allow cross-origin requests.

## Supported files

- Mantid `MDHistoWorkspace` (`/MDHistoWorkspace/data/signal`, `mask`, `D0`–`D2`). The cell comes from the UB matrix in `experiment0/sample/oriented_lattice/orientation_matrix` (reciprocal metric = UBᵀ·UB), falling back to `unit_cell_*`. Axis names such as `[H,0,0]` or `[H,H,0]` give the basis vectors used for the oblique geometry.
- Any other `NXdata` group with a 3-D signal (found through `@signal`, and `@axes` on the group or the signal). Axes may hold bin edges or bin centers. Size-1 dimensions are ignored, so a 4-D workspace with one integrated axis works. Without an HKL frame and a unit cell, axes are drawn rectangular.
- Compression: deflate (gzip), shuffle, and the other filters built into HDF5. Plugin filters (LZ4, Blosc, bitshuffle) are not supported yet.

## Symmetry operations

Operations are written like coordinate triplets and act on reflections: `h+k,-h,l` is the six-fold rotation about c* in hexagonal axes, `-h,-k,-l` is inversion. Real-space triplets such as `x-y,x,z` from the International Tables are accepted and converted (a reflection h is equivalent to Wᵀh). Separate operations with `;`. The viewer closes the generators into a group (at most 48 operations) and lists every operation.

Presets: 1, −1, 2/m (b or c unique), mmm, 4/m, 4/mmm, −3, −3m1, −31m, 6/m, 6/mmm, m−3, m−3m (trigonal and hexagonal classes use hexagonal axes).

Two checks guard against mistakes:

- **Grid**: each operation must send bin centers onto bin centers, so equivalent voxels are found exactly, without interpolation. Coupled axes need equal bin widths, and rotations other than inversion need a bin centered at the origin. Operations use the HKL basis of each axis (for example `[H,H,0]`), so they also work on projected grids.
- **Metric**: the viewer reports how much the operations change the reciprocal metric of the cell. It warns above 2%, which usually means the operations belong to a different setting. Fe₃Ge₂ is pseudo-hexagonal, and 6/mmm changes its measured metric by 0.8%.

## How slices are computed

Each output pixel is the equal-weight mean of the unique finite, unmasked source voxels in the slab and all their symmetry equivalents. Symmetry orbits are either identical or disjoint, so the pooled set is the union of the distinct orbits met in the slab column, and each voxel is counted once, including voxels on symmetry elements. Voxels with a nonzero `mask` or a non-finite signal are excluded. Grey pixels have no measured contribution. There is no smoothing or interpolation.

With no symmetry or with −1 this is exactly `average_slab()` from the original Fe₃Ge₂ viewer. The port is checked against that function in the test suite, and on a real 401³ file counts and empty pixels matched exactly and values agreed to within float32 precision. Every Laue class is also checked against a direct set-based implementation. A 3-bin slice with 6/mmm takes about 40 ms.

The 3-D view bins the volume by an odd factor, with blocks aligned to the origin so the same operations map blocks onto blocks. It pools block sums and counts over each orbit and draws the isosurface of the block means with surface nets. Its default level is the 99.5th percentile of the positive block means.

The signal is shown as stored in the file (no `num_events` normalization). Voxels are held as float32 and sums are accumulated in float64. Non-uniform bins are drawn as if uniform.

## Development

There is no build step. Serve the folder with any static server and open it:

```bash
python3 -m http.server 8000
```

Tests run in Node against small fixture files:

```bash
npm install && npm test
```

- `js/slab.js`: slab selection and symmetry-pooled averaging (pure functions).
- `js/symmetry.js`: parsing operations, group closure, metric check, and index maps on the bin grid.
- `js/iso.js`: coarse binning, orbit means, and surface nets for the 3-D view.
- `js/nexus.js`: finds the histogram, reads it with h5wasm, and handles UB and lattice geometry.
- `js/worker.js`: module worker that holds the volume and answers slice and isosurface requests.
- `js/view3d.js`: three.js scene, loaded when the 3-D panel opens.
- `js/app.js`: UI and canvas rendering.
- `tools/make_fixtures.py`: regenerates `tests/fixtures` and expected values from the reference `average_slab()` (needs numpy and h5py).
- `tools/make_colormaps.py`: regenerates `js/colormaps.js` from matplotlib.

h5wasm and three.js are loaded from jsDelivr at pinned versions (`js/worker.js` and the import map in `index.html`); keep h5wasm in step with `package.json`.

## Publishing on GitHub Pages

In the repository settings, go to **Pages** and choose **Deploy from a branch**, then select `main` and `/ (root)`. The site will be at `https://<user>.github.io/neutron-nexus-viewer/`.

## Credits

HDF5 reading uses [h5wasm](https://github.com/usnistgov/h5wasm) from the National Institute of Standards and Technology. 3-D rendering uses [three.js](https://threejs.org).
