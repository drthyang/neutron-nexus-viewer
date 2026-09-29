# User guide

[← README](../README.md) · [Method](METHOD.md) · [Development](DEVELOPMENT.md)

## Opening data

- **Local file**: click the dataset button or drop a `.nxs` file anywhere on the page. The file is read in the browser and never uploaded.
- **Example**: *Try the example* on the start page (or `?demo`) opens one synthetic hexagonal crystal at two temperatures, compared in split view with 6/mmm averaging and a 1-voxel edge mask. At 300 K (`examples/demo_300K.nxs`), short-range order gives diffuse rods along L at the M points; at 10 K (`examples/demo_10K.nxs`), they condense into superlattice peaks at even L. Both have coverage gaps, which symmetry averaging fills, and bright detector-edge voxels, which the mask removes: set the Laue class to *None* and *Clear* the mask to see the raw data.
- **Links**: `?url=https://…/file.nxs` downloads and opens a file, and `&compare=https://…/other.nxs` adds a second file to compare. `&sym=6/mmm` applies a Laue class, and `&mask=1` (or `&mask=1,5`) applies a mask with that erosion radius (and outlier cut). The host must allow cross-origin requests.

The histogram is held in memory as float32. A 401³ volume needs about 260 MB and takes 1–2 s to read.

## Supported files

- **Mantid `MDHistoWorkspace`** (`SaveMD`): `/MDHistoWorkspace/data/signal`, `mask` and `D0`–`D2`. The cell comes from the UB matrix in `experiment0/sample/oriented_lattice/orientation_matrix`, falling back to `unit_cell_*`. Axis names such as `[H,0,0]` or `[H,H,0]` give the basis vectors for the oblique geometry.
- **Any other `NXdata` group with a 3-D signal**, found through `@signal`, and `@axes` on the group or the signal. Axes may hold bin edges or bin centers. Size-1 dimensions are ignored, so a 4-D workspace with one integrated axis works. Without an HKL frame and a unit cell, the axes are drawn rectangular.
- **Compression**: deflate (gzip), shuffle and the other filters built into HDF5. Plugin filters (LZ4, Blosc, bitshuffle) are not supported.
- **Limits**: 3-D histograms up to about 3.5 GB as float32. The signal is shown as stored, with no `num_events` normalization. Non-uniform bins are drawn as if uniform.

## The screen

- **Top bar**: the open file (A and B, joined by the split icon, when comparing), the engine status, and ⓘ for dataset details, a summary of the method and links to these guides. Long file names are shortened in the middle; hover for the full name.
- **Control panel** (left, collapsible; each section folds to a one-line summary), in the order you work:
  - **Dataset**: unit cell (from UB), reciprocal lattice, grid and measured fraction, and the Compare card for opening a second dataset.
  - **Display**: colormap, color range, scale, view range, guides and cell angles.
  - **Processing**: the pipeline *Measured voxels → Mask → Symmetry average → Views → NEBULA3D*, followed by the Symmetry and Mask controls. Active stages are highlighted, and clicking a stage jumps to its controls.
  - **Export**: the last step, handing the processed volume to NEBULA3D.
- **Workspace**: a header row with the click mode, the shared color legend (click it for the Display settings) and the layout, then four views (HK, HL, KL and 3-D).
- **Panel footer**: copyright, the license and a link to this documentation. The ⓘ popover links to each guide.

Below 1000 px wide, or on short screens, the panel sits above the views and the page scrolls.

## Symmetry averaging

Pick a Laue class, or type generators and press *Apply*. Operations are written like coordinate triplets and act on reflections: `h+k,-h,l` is the six-fold rotation about c* in hexagonal axes, and `-h,-k,-l` is inversion. Real-space triplets such as `x-y,x,z` are also accepted. Separate operations with `;`. *Group operations* lists the closed group.

The status line reports how much the operations change the cell metric. A warning above 2% usually means the operations belong to a different setting. Operations that do not map bin centers onto bin centers are rejected with an explanation. See [Method → Symmetry](METHOD.md#symmetry-averaging).

## Mask

The mask removes voxels from the measured data before symmetry averaging:

- **Edge erosion**: removes measured voxels within *r* voxels of unmeasured ones. Detector edges leave spurious values there.
- **Outlier cut k·σ**: removes voxels more than *k* robust standard deviations above the median of their symmetry equivalents. It needs a Laue class; 0 disables it.

*Apply mask* builds it (a few seconds for a 401³ volume), and *Clear* removes it. *Show removed voxels* averages only what the mask took out, so you can check that it hits artifacts. *Download* saves the mask as a gzipped NumPy `.npy` in the signal's storage order (1 = edge, 2 = outlier): `signal[mask > 0] = np.nan` applies it in Python. See [Method → Mask](METHOD.md#masking-detector-edge-artifacts).

## Display

- **Colormap and range**: vmin and vmax set the ends of the colormap. The legend above the views shows the current scale.
- **Scale**: asinh (with a softening value), linear or log. *Auto range* sets vmin to 0, vmax to the 97th percentile, and softening to the median of the positive values in the current slices.
- **View ±**: shows ± this many r.l.u. around the origin. Leave it empty for the full range.
- **Guides**: dashed lines where the other two slices cut each view.
- **Nominal cell angles**: snaps direct-cell angles within 1° of 60°, 90° or 120° for drawing. Otherwise the angles derived from UB are used.

Axes are drawn with the reciprocal metric of the cell. Axes that share a unit have the same length per r.l.u. and share one tick step. Views whose axes are not orthogonal show the angle between them, for example ∠ 60° for HK in a hexagonal cell.

## Slice views

- **Footer**: slab center and full thickness of the integrated axis, as slider or number. Bins whose centers lie within thickness/2 of the center are averaged. Views update live while you drag.
- **Header**: plane, position, coverage and processing, with buttons for *Reset zoom*, *Save PNG*, *focus* and *maximize*.
- **Hover**: shows the coordinates, the value and the number of pooled voxels.

### Clicking a slice

The workspace header switches what the mouse does:

| Mode | Click | Drag |
| --- | --- | --- |
| **Navigate** | moves the other two slices through the point | the same, continuously |
| **Zoom** | zooms in 2× around the point | zooms into the box; with a common unit, the window is widened so both axes span the same length |
| **Move** | — | slides the visible region |

In every mode, double-click (or *Reset zoom*) returns to the full view.

## 3-D view

A transparent isosurface of the binned, symmetrized and masked volume, with the current slices as planes, clipped to the view range. The footer sets the isosurface level (log slider or typed; empty returns to the automatic level) and the surface and slice opacity. The options button sets the grid (about 64, 100 or 150 blocks per axis) and hides the slices. Drag to rotate, scroll to zoom, right-drag to pan, and use ↺ to reset the camera.

## Comparing two datasets

*Open second file (B)* in the Compare card, or drop a file on the card, opens a second dataset next to the first (A). Every slice is then cut along its diagonal, from the top-left to the bottom-right corner of the view: A fills the lower-left half and B the upper-right half. A white gap marks the cut, B's half is hatched where it has no data, and tags in the corners name the files.

Once B is open, the Dataset section lists both files and one table for the two: a value they share (often the cell and grid) appears once, and a value that differs gets a line for each, tagged A and B.

- **Shared**: slice positions and thickness, zoom and pan, symmetry, mask settings and the color scale apply to both datasets, so the two halves are directly comparable. *Auto range* pools both.
- **A / Split / B** in the workspace header shows one dataset over the whole view, or the split.
- **Hover** reads both datasets at the cursor, the one under it first. View headers show the coverage of each.
- **Processing**: B has its own worker and mask, built with the same parameters on its own data. Symmetry operations are mapped onto B's grid; if they do not fit it, B is used as measured and the card says so.
- **Axes**: B is drawn on A's axes and lattice geometry. Files should share the same axes (for example both `[H,0,0]`, `[0,K,0]`, `[0,0,L]`); the grids may differ. A warning under the Dataset table appears when the axis names differ.
- **3-D view**: the isosurface is A's; the slice planes show the same split as the views.
- **Replace** (⇄ on B's line, or the B half of the top bar) swaps in another file, and **Remove** (×) returns to a single dataset. Opening a new file as A also removes B.

## Export for NEBULA3D

The **Export** section at the bottom of the panel hands the processed volume to the 3D-ΔPDF pipeline of [NEBULA3D](https://github.com/drthyang/nebula3d). NEBULA3D does not symmetrize, so choose a Laue class (and a mask) first.

- **Open in NEBULA3D** opens NEBULA3D in a new tab and builds the volume while it starts; NEBULA3D then loads it and selects it as its dataset, ready to configure and run. The card reports when NEBULA3D has loaded it. If the browser blocks the tab, allow pop-ups for this site.
- **Download** saves the same file, to open in NEBULA3D with *Load volume…* (or to put in the data folder of its desktop app).

- **Content**: every voxel is the mean of the measured, unmasked voxels among its symmetry equivalents, as in the slices; voxels with none are written as 0 with mask 0, and NEBULA3D backfills them. When comparing, dataset A is exported.
- **Grid**: axes must be plain H, K and L with uniform bins. Each axis is padded to be symmetric about 0 (NEBULA3D puts Q = 0 at the centre), and symmetry fills the padding where equivalents were measured.
- **File**: HDF5 with `/entry/{data, mask, h_axis, k_axis, l_axis, ub_matrix}` in (H, K, L) order, named like `<file>_sym6mmm.nxs`, with the source file, symmetry, mask and UB source as attributes. NEBULA3D's browser build handles up to about 80 M voxels; a 401³ volume (64 M) takes about 10 s to export. See [Method → Export](METHOD.md#export-for-nebula3d).

## Layouts and export

- **Layouts**: *Quad* (2×2), *Focus* (one large view with the other three beside it) and *Single*. In *Focus*, the small views are thumbnails: hovering one highlights it, and clicking it (or Enter) shows it large. Each view's header can focus or maximize it, and double-clicking a header maximizes it. Esc returns.
- **Save PNG** exports a view at 3× resolution. Slice exports include a title and their own colorbar, without guides. When comparing, they keep the split and the dataset tags.
- **Remembered settings**: the colormap, scale, click mode, layout, panel state and folded sections are remembered per browser.
