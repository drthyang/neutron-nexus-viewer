# NeXus slice viewer

**Live site: [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/)**

A static web page for browsing reciprocal-space slices of 3-D neutron histograms, such as Mantid `MDHistoWorkspace` files written by `SaveMD`. Open or drop a `.nxs` file and the page shows three orthogonal slices (HK, HL, KL) with adjustable slab center and thickness, a 3-D view with a transparent isosurface and the current slices, symmetry averaging with any Laue class or your own operations, and masking of detector-edge artifacts before averaging. The interface follows the nebula3d console design, in a light theme, on a single page.

The file is read in your browser by [h5wasm](https://github.com/usnistgov/h5wasm) (HDF5 compiled to WebAssembly) and is never uploaded. The page is static, so it can be hosted on GitHub Pages.

## Using it

The page is a full-window workspace: a top bar (dataset, symmetry and mask chips, layout, info), a toolbar with the shared color scale, and four views (HK, HL, KL and 3-D). Below 1000 px wide, or on short screens, the views stack and the page scrolls.

- **Open a file** with the dataset button or drag it onto the page. The histogram is loaded into memory as float32: a 401³ volume needs about 260 MB and takes 1–2 s to read.
- **Layout**: *Quad* shows the four views in a 2×2 grid; *Focus* enlarges one view with the other three stacked beside it; *Single* shows one view. Each view's header has buttons to focus or maximize it; double-clicking a header maximizes it, and Esc returns.
- **Slices**: each view's footer sets the slab center and full thickness of its integrated axis. Bins whose centers lie within thickness/2 of the center are averaged, exactly as in the original Python viewer. Slices update live while you drag. Hovering shows the coordinates, value and pooled voxel count.
- **Clicking a slice** has two modes, switched in the toolbar:
  - *Navigate* (default): click, or drag, to move the other two slices so all three intersect at the pointer. Dashed guides show where they cut.
  - *Zoom*: drag a rectangle to zoom into it (with oblique axes the window is the H–K range enclosing the box, previewed while dragging); a plain click zooms in 2× around the point.
  - In both modes, double-click or *Reset zoom* returns to the full view.
- **Color** (toolbar): colormap, vmin/vmax around the shared colorbar, asinh, linear or log scale, and asinh softening. *Auto range* sets vmax to the 97th percentile and softening to the median of the positive values in the current slices.
- **View** (toolbar): view range ±, dashed guides, and *Nominal cell angles*. Axes are drawn with the reciprocal metric from the UB matrix in the file (or its stored unit cell when there is no UB); nominal angles snap direct-cell angles within 1° of 60°, 90° or 120° (for example, 90/90/120 for a hexagonal cell).
- **Symmetry** (top-bar chip): pick a Laue class or type generators (see below). Slices and the 3-D view pool every voxel with its symmetry equivalents. The default is no symmetry, so the data are shown as measured.
- **Mask** (top-bar chip): removes spurious voxels from the unsymmetrized data before averaging (see below). *Show removed* averages only the removed voxels, so you can check what the mask takes out. *Download* saves the mask as a gzipped NumPy `.npy`.
- **3-D view**: a transparent isosurface of the binned volume with the three current slices as planes, drawn in the lattice geometry and clipped to the view range. The footer sets the level (log slider or typed) and the surface and slice opacity (lower slice opacity makes the planes see-through); the options button sets the grid (about 64, 100 or 150 blocks per axis) and hides the slices. Drag to rotate, scroll to zoom, right-drag to pan.
- **Info** (ⓘ): file, axes, grid, cell, measured fraction, current symmetry and mask, and a summary of the method.
- **Save PNG** (camera button) exports a view at 3× resolution; slice exports include a title and their own colorbar, without guides.
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

## Masking detector-edge artifacts

Detector edges, and the weak normalization there, leave spuriously high values along the boundaries of the measured region in reciprocal space. The reduced file no longer knows which detector pixel a voxel came from. Two masks act on the unsymmetrized volume instead, and symmetry averaging then uses only the voxels that remain:

- **Edge erosion (voxels)**: removes measured voxels within *r* voxels (box distance) of an unmeasured voxel. In the Fe₃Ge₂ 90 K volume, the 99th percentile of the voxels within one voxel of an edge is 4266, against 36.5 for interior voxels. *r* = 1 removes 15.8% of the measured voxels and *r* = 2 removes 30.6%.
- **Outlier cut (k·σ)**: for each symmetry orbit with at least 3 valid voxels, removes voxels more than *k* robust standard deviations (1.4826 × MAD) above the orbit median. This needs a Laue class to be selected. With 6/mmm and *k* = 5 it removes about 1% of voxels.

With *r* = 2 and *k* = 5, the 6/mmm-averaged HK plane at L = 0 has 120 pixels above 1000 instead of 1476, and its coverage only drops from 73.9% to 71.7%, because symmetry fills most gaps. The 99.5th percentile of the 3-D block means falls from 4340 to 88: most of the brightest blocks were edge artifacts.

The downloaded mask is `uint8` in the signal dataset's storage order (1 = edge, 2 = outlier), so `signal[mask > 0] = np.nan` applies it in Python. The cleaner fix is upstream: mask detector-edge pixels (for example with Mantid's `MaskBTP`) before converting to MD and normalizing.

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
- `js/mask.js`: coverage-edge erosion and symmetry-outlier masks.
- `js/nexus.js`: finds the histogram, reads it with h5wasm, and handles UB and lattice geometry.
- `js/worker.js`: module worker that holds the volume and answers slice and isosurface requests.
- `js/view3d.js`: three.js scene, loaded when the 3-D panel opens.
- `js/app.js`: UI and canvas rendering.
- `tools/make_fixtures.py`: regenerates `tests/fixtures` and expected values from the reference `average_slab()` (needs numpy and h5py).
- `tools/make_colormaps.py`: regenerates `js/colormaps.js` from matplotlib.

h5wasm and three.js are loaded from jsDelivr at pinned versions (`js/worker.js` and the import map in `index.html`); keep h5wasm in step with `package.json`.

## Publishing on GitHub Pages

The site is published by GitHub Pages from the root of `main` (**Deploy from a branch**, `main`, `/ (root)`) at [drthyang.github.io/neutron-nexus-viewer](https://drthyang.github.io/neutron-nexus-viewer/). Every push to `main` redeploys it, usually within a minute. A fork can publish its own copy the same way, at `https://<user>.github.io/neutron-nexus-viewer/`.

## Credits

HDF5 reading uses [h5wasm](https://github.com/usnistgov/h5wasm) from the National Institute of Standards and Technology. 3-D rendering uses [three.js](https://threejs.org).
