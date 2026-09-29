# User guide

[← README](../README.md) · [Method](METHOD.md) · [Development](DEVELOPMENT.md)

## Opening data

- **Local file**: click the dataset button or drop a `.nxs` file anywhere on the page. The file is read in the browser and never uploaded.
- **Example**: *Try the example* on the start page (or `?demo`) opens one synthetic hexagonal crystal at two temperatures, compared in split view with 6/mmm averaging and a 1-voxel edge mask. At 300 K (`examples/demo_300K.nxs`), short-range order gives diffuse rods along L at the M points; at 10 K (`examples/demo_10K.nxs`), they condense into superlattice peaks at even L. Both have coverage gaps, which symmetry averaging fills, and bright detector-edge voxels, which the mask removes: set the Laue class to *None* and *Clear* the mask to see the raw data.
- **Links**: `?url=https://…/file.nxs` downloads and opens a file, and `&compare=https://…/other.nxs` adds a second file to compare. `&sym=6/mmm` applies a Laue class, and `&mask=1` (or `&mask=1,5`) applies a mask with that erosion radius (and outlier cut). The host must allow cross-origin requests.

The histogram is held in memory as float32. A 401³ volume needs about 260 MB and takes 1–2 s to read.

## Supported files

- **Mantid `MDHistoWorkspace`** (`SaveMD`): `/MDHistoWorkspace/data/signal`, `mask` and `D0`–`D2`, and `errors_squared` for the uncertainties of I(Q). The cell comes from the UB matrix in `experiment0/sample/oriented_lattice/orientation_matrix`, falling back to `unit_cell_*`. Axis names such as `[H,0,0]` or `[H,H,0]` give the basis vectors for the oblique geometry.
- **Any other `NXdata` group with a 3-D signal**, found through `@signal`, and `@axes` on the group or the signal. Axes may hold bin edges or bin centers. Size-1 dimensions are ignored, so a 4-D workspace with one integrated axis works. Without an HKL frame and a unit cell, the axes are drawn rectangular.
- **Compression**: deflate (gzip), shuffle and the other filters built into HDF5. Plugin filters (LZ4, Blosc, bitshuffle) are not supported.
- **Limits**: 3-D histograms up to about 3.5 GB as float32. The signal is shown as stored, with no `num_events` normalization. Non-uniform bins are drawn as if uniform.

## The screen

- **Top bar**: the datasets, then the status and ⓘ for the full dataset details, a summary of the method and links to these guides. The status also shows work in progress (building a mask, I(Q), the NEBULA3D export, opening B): the job that started first and how far it is, with the number of others ("Mask · 45% +1"), and a thin bar along the bottom edge of the top bar for all of them together; hover it for the list. Each dataset is a chip with its file name over one line of facts: the cell (angles other than 90°), the grid and the measured fraction; hover for the details, such as the reciprocal lattice, bin widths and mask. Click A to open another file, or *Compare…* to add a second one (B); see [Comparing](#comparing-two-datasets). Long file names are shortened in the middle.
- **Control panel** (left, collapsible; each section folds to a one-line summary), in the order you work:
  - **Processing**: the pipeline *Measured voxels → Mask → Symmetry average → Views → NEBULA3D*, followed by the Symmetry and Mask controls. Active stages are highlighted, and clicking a stage jumps to its controls.
  - **Export**: the last step: I(Q), the volume reduced to 1-D, and handing the processed volume to NEBULA3D.
- **Workspace**: a header row with the click mode, the display controls (see [Display](#display)) and the layout, then four views: HK, HL, KL, and the 3-D view or I(Q), which share the fourth place (the *3D | I(Q)* switch in its header picks one).
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

The color scale is set in the workspace header, above the views: the colormap, vmin and vmax on either side of the colorbar, and the scale.

- **Scale**: asinh (with a softening value), lin or log. *Auto* sets vmin to 0, vmax to the 97th percentile, and softening to the median of the positive values in the current slices.

The sliders button next to *Auto* opens the other options:

- **Softening** (asinh only): values below it are shown nearly linearly.
- **View ±**: shows ± this many r.l.u. around the origin. Leave it empty for the full range.
- **Guides**: dashed lines where the other two slices cut each view.
- **Integer grid**: thin lines at whole-number values of each view's axes, such as integer H and K in the HK view. They follow the true axis directions, so on a hexagonal cell they meet at 60°. An axis whose lines would be closer than 6 pixels gets none.
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

## I(Q)

I(Q) reduces the masked, symmetrized volume to one dimension: the mean intensity in each shell of |Q|, over the part of the shell that has data. Open it with *I(Q)* in the header of the 3-D view (the two share the fourth place) or *Show I(Q)* in the **Export** section. It is computed when it is shown, and again when the symmetry, the mask or its settings change.

- **Normalization**: unmeasured and masked voxels are left out, not counted as zero, so gaps in coverage do not lower I(Q). Symmetry-equivalent voxels are pooled as in the slices, and each orbit counts with its multiplicity, so a shell is not biased toward the directions that were measured best. See [Method → I(Q)](METHOD.md#powder-average-iq).
- **Q bins** (footer): the ΔQ field takes Mantid `Rebin` parameters. `0.05` gives shells of 0.05 Å⁻¹; a negative step gives logarithmic shells, so `-0.01` means ΔQ/Q = 1% (suited to the constant Δd/d resolution of time-of-flight instruments); `0.5, 0.02, 3, 0.05, 10` gives 0.02 Å⁻¹ shells from 0.5 to 3 Å⁻¹ and 0.05 Å⁻¹ shells from 3 to 10 Å⁻¹, and any mix of ranges works. As in Mantid, a range ends at its boundary with a last bin of 0.25 to 1.25 steps. Empty means shells of the shortest bin step in |Q|.
- **ΔQ slider**: changes a single step on a log scale, from a tenth of the shortest bin step to 20 times it (or ΔQ/Q from 0.1% to 20% for logarithmic shells), and recomputes I(Q) as it moves. Typed values apply on Enter.
- **Q range**: Q min (empty for 0; the shortest bin step for logarithmic shells) and Q max (empty for all the data) bound a single step; ranges typed in the ΔQ field set their own, and the fields are then disabled. *Lin | Log* in the view's header switches the intensity scale. |Q| is in Å⁻¹ with 2π, from the file's UB matrix or cell; axes already in Å⁻¹ (Mantid's Q frames) need no cell.
- **Options** (the sliders button): *Split voxels* shares each voxel between the shells it overlaps by dividing it into 2³ or 3³ sub-cells, each binned by its own |Q|; *Centres* bins whole voxels, which is faster but aliases when ΔQ is close to the voxel size. The error band (±σ) and the shell coverage (the dashed line, right axis) can be hidden.
- **Uncertainties**: σ is propagated from the file's `errors_squared` (or `errors`), read the first time I(Q) is computed. Without them, the curve has no error band and the text file has `nan` in its σ column.
- **Coverage**: the fraction of each shell's volume that has data, after symmetry. Where it falls (beyond the measured region, or where shells leave the grid), I(Q) rests on few voxels.
- **Plot**: hover to read Q, d = 2π/Q, and I ± σ with the coverage of each dataset. Drag a box to zoom (a flat drag zooms Q only), click in *Zoom* mode to zoom 2×, drag in *Move* mode to pan, and double-click (or *Reset zoom*) for the full range.
- **Comparing**: both datasets are reduced on the same shells, each with its own grid, cell and mask, and drawn as two curves; *A / Split / B* chooses which are shown.
- **Download** (in the view's header or the Export section) saves a text file, `<file>_IQ.dat`: a commented header with the file, cell, shells (also as Rebin parameters, to reproduce them in Mantid), symmetry, mask and normalization, then one row per shell with Q (the shell centre, the midpoint of its edges), I, σ, coverage and the number of voxels with data, for A and B when comparing. `numpy.loadtxt` reads it directly. *Save PNG* exports the plot at 3× with a title and a legend.

## Comparing two datasets

*Compare…* in the top bar, or dropping a file on it, opens a second dataset next to the first (A). Every slice is then cut along its diagonal, from the top-left to the bottom-right corner of the view: A fills the lower-left half and B the upper-right half. A white gap marks the cut, B's half is hatched where it has no data, and tags in the corners name the files.

B then gets its own chip in the top bar, after A and the split icon, with its facts. While it loads, and while its mask is built, the chip shows the step, and the status in the top bar the progress; an amber *!* marks warnings (hover to read them).

- **Shared**: slice positions and thickness, zoom and pan, symmetry, mask settings and the color scale apply to both datasets, so the two halves are directly comparable. *Auto range* pools both.
- **A / Split / B** in the workspace header shows one dataset over the whole view, or the split.
- **Hover** reads both datasets at the cursor, the one under it first. View headers show the coverage of each.
- **Processing**: B has its own worker and mask, built with the same parameters on its own data. Symmetry operations are mapped onto B's grid; if they do not fit it, B is used as measured and the card says so.
- **Axes**: B is drawn on A's axes and lattice geometry. Files should share the same axes (for example both `[H,0,0]`, `[0,K,0]`, `[0,0,L]`); the grids may differ. B's chip warns when the axis names differ.
- **3-D view**: the isosurface is A's; the slice planes show the same split as the views.
- **Replace**: click B's chip to swap in another file. **Remove** (× on the chip, which also cancels a loading B) returns to a single dataset. Opening a new file as A also removes B.

## Export for NEBULA3D

The **Export** section at the bottom of the panel hands the processed volume to the 3D-ΔPDF pipeline of [NEBULA3D](https://github.com/drthyang/nebula3d). NEBULA3D does not symmetrize, so choose a Laue class (and a mask) first.

- **Open in NEBULA3D** opens NEBULA3D in a new tab and builds the volume while it starts, and NEBULA3D shows the build's progress meanwhile; it then loads the volume and selects it as its dataset, ready to configure and run. The card reports when NEBULA3D has loaded it. If the browser blocks the tab, allow pop-ups for this site.
- **Download** saves the same file, to open in NEBULA3D with *Load volume…* (or to put in the data folder of its desktop app).

- **Content**: every voxel is the mean of the measured, unmasked voxels among its symmetry equivalents, as in the slices; voxels with none are written as 0 with mask 0, and NEBULA3D backfills them. When comparing, dataset A is exported.
- **Grid**: axes must be plain H, K and L with uniform bins. Each axis is padded to be symmetric about 0 (NEBULA3D puts Q = 0 at the centre), and symmetry fills the padding where equivalents were measured.
- **File**: HDF5 with `/entry/{data, mask, h_axis, k_axis, l_axis, ub_matrix}` in (H, K, L) order, named like `<file>_sym6mmm.nxs`, with the source file, symmetry, mask and UB source as attributes. NEBULA3D's browser build handles up to about 80 M voxels; a 401³ volume (64 M) takes about 10 s to export. See [Method → Export](METHOD.md#export-for-nebula3d).

## Layouts and export

- **Layouts**: *Quad* (2×2), *Focus* (one large view with the other three beside it) and *Single*. In *Focus*, the small views are thumbnails: hovering one highlights it, and clicking it (or Enter) shows it large. Each view's header can focus or maximize it, and double-clicking a header maximizes it. Esc returns.
- **Save PNG** exports a view at 3× resolution. Slice exports include a title, their own colorbar and the integer grid when it is on, without guides. When comparing, they keep the split and the dataset tags.
- **Remembered settings**: the colormap, scale, click mode, layout, the view in the fourth place, the I(Q) options, panel state and folded sections are remembered per browser.
