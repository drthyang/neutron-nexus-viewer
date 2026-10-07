# User guide

[← README](../README.md) · [Method](METHOD.md) · [Development](DEVELOPMENT.md)

## Opening data

- **Local file**: click the dataset button or drop a `.nxs` or `.h5` file anywhere on the page. The file is read in the browser and never uploaded.
- **Example**: *Try the example* on the start page (or `?demo`) opens one synthetic hexagonal crystal at two temperatures, compared in split view with 6/mmm averaging and a 1-voxel edge mask. At 300 K (`examples/demo_300K.nxs`), short-range order gives diffuse rods along L at the M points; at 10 K (`examples/demo_10K.nxs`), they condense into superlattice peaks at even L. Both have coverage gaps, which symmetry averaging fills, and bright detector-edge voxels, which the mask removes: set the Laue class to *None* and *Clear* the mask to see the raw data.
- **Rigaku XRD frames**: *Reduce Rigaku XRD…* on the start page, or the reduction button in the top bar, turns a CrysAlisPro experiment folder into an HKL volume in the browser; see [Reducing Rigaku XRD frames](#reducing-rigaku-xrd-frames).
- **Links**: `?url=https://…/file.nxs` downloads and opens a file, and `&compare=https://…/other.nxs` adds a file to compare (repeat it, up to three times, for B, C and D). `&sym=6/mmm` applies a Laue class, and `&mask=1` (or `&mask=1,5`) applies a mask with that erosion radius (and outlier cut). The host must allow cross-origin requests.

The histogram is held in memory as float32. A 401³ volume needs about 260 MB and takes 1–2 s to read.

## Supported files

- **Mantid `MDHistoWorkspace`** (`SaveMD`): `/MDHistoWorkspace/data/signal`, `mask` and `D0`–`D2`, and `errors_squared` for the uncertainties of I(Q). The cell comes from the UB matrix in `experiment0/sample/oriented_lattice/orientation_matrix`, falling back to `unit_cell_*`. Axis names such as `[H,0,0]` or `[H,H,0]` give the basis vectors for the oblique geometry.
- **Any other `NXdata` group with a 3-D signal**, found through `@signal`, and `@axes` on the group or the signal. Axes may hold bin edges or bin centers. Size-1 dimensions are ignored, so a 4-D workspace with one integrated axis works. Without an HKL frame and a unit cell, the axes are drawn rectangular. Axes named `x`, `y` and `z` in Å are taken as real space along the cell axes a, b and c (a 3D-ΔPDF; see below).
- **[NEBULA3D](https://github.com/drthyang/nebula3d) files** written before it adopted the Mantid layout, which have no `NXdata`:
  - *Volumes* (`_ringremoved.h5`, `_braggpunched.h5`, `_backfilled.h5`, `_flattened.h5`): `/entry/data` in (H, K, L) order, `mask` (1 = valid), `sigma` (standard deviations, squared for I(Q) and line cuts), bin-center axes `h_axis`, `k_axis`, `l_axis`, and `ub_matrix`, which includes 2π.
  - *3D-ΔPDFs* (`_delta_pdf.h5`): `/data` with `x_axis`, `y_axis`, `z_axis` in Å along a, b and c, and the cell in the `lat_a` … `lat_gamma` attributes. Sections are drawn at the cell angles, and symmetry operations act on x, y, z as the real-space form of the given h, k, l operations. A ΔPDF changes sign, so it opens with the coolwarm colormap on a color range symmetric about 0; RdBu_r, also diverging, is in the colormap list. I(Q) and the NEBULA3D export need reciprocal-space axes and are not available for it.

  As with Mantid's D0, D1, D2, the axes are listed fastest-varying first: L, K, H for volumes, z, y, x for ΔPDFs.
- **Compression**: deflate (gzip), shuffle and the other filters built into HDF5. Plugin filters (LZ4, Blosc, bitshuffle) are not supported.
- **Limits**: 3-D histograms up to about 3.5 GB as float32. The signal is shown as stored, with no `num_events` normalization. Non-uniform bins are drawn as if uniform.

## Reducing Rigaku XRD frames

*Reduce Rigaku XRD…* (on the start page, or the button next to ⓘ in the top bar) makes an HKL volume from the raw frames of a Rigaku Oxford Diffraction single-crystal experiment, as CrysAlisPro writes them.

1. **Choose folder…** and pick the experiment folder: the one holding `frames/` with the `*.rod_img` files, the `.par` files and `expinfo/`. The browser asks to "upload" the folder; nothing leaves your computer. The dialog lists:
   - the runs, each with its scan, κ, φ and exposure;
   - the detector, wavelength, monochromator, temperature and Laue class;
   - the cells of the orientation matrices CrysAlisPro stored;
   - whether the first frame decodes to the statistics in its header.
2. **Choose the options:**
   - **Runs**: untick runs to leave them out.
   - **Output cell**: the indexing of the volume. *CrysAlis cell* indexes in the refined cell. *2 × 2 × 2 cell* doubles every index, as in neutron reductions that use a doubled cell. *Match dataset A* (when a volume is open) takes the multiple from A's cell and puts the voxel centres on A's, so the two compare voxel by voxel.
   - **Voxel**: the bin width in output r.l.u. The default, 0.05 in the CrysAlis cell, is about the width of the Bragg peaks and of a 0.5° frame at high Q.
   - **Steps per frame**: each frame's rotation is split into this many steps sharing its counts.
   - **Normalization**: *Solid angle + polarization* (the default) divides by each pixel's solid angle and the polarization factor. *Exposure only* gives counts per second per pixel.
   - **Refine the geometry**: fits the detector, goniometer offsets and orientation (per run) to the Bragg peaks in the frames, starting from the best CrysAlisPro matrix. Turn it off only to see the stored model as it is.
   - **Open the result**: as a new dataset, or to compare with the open one.
3. **Reduce.** The frames are read three times: for the detector mask, for the Bragg peaks, and for the mapping. Expect about half a minute per thousand frames on a recent laptop (3,608 frames took 2 minutes) and 1–1.5 GB of memory. *Cancel* stops it.

The report lists:
- the frame checks and the mask;
- the peaks found and indexed, and with which matrix;
- the refined geometry's residuals and any runs that needed a drift model;
- the cell, the grid and the measured voxels.

*Download .nxs* saves a Mantid `MDHistoWorkspace` that `LoadMD` reads. Voxels without data are NaN, and `num_events` counts the detector pixel-frames in each voxel. *Report (.json)* saves the full report with the refined model.

Not applied: background subtraction (air scatter and fluorescence show near Q = 0), absorption (it needs the crystal's shape), and symmetry averaging (use **Symmetry averaging** in the panel). Where runs with different absorption meet, the volume can show steps. The method is in [Method → Rigaku reduction](METHOD.md#rigaku-reduction).

## The screen

- **Top bar**: the datasets, then the status and ⓘ for the full dataset details, a summary of the method and links to these guides. The status also shows work in progress (building a mask, I(Q), the NEBULA3D export, opening B, C or D): the job that started first and how far it is, with the number of others ("Mask · 45% +1"), and a thin bar along the bottom edge of the top bar for all of them together; hover it for the list. Each dataset is a chip with its file name over one line of facts: the cell (angles other than 90°), the grid and the measured fraction; hover for the details, such as the reciprocal lattice, bin widths and mask. Click A to open another file, or *Compare…* to add another one (B, C, D); see [Comparing](#comparing-datasets). Long file names are shortened in the middle.
- **Control panel** (left, collapsible; each section folds to a one-line summary), in the order you work:
  - **Processing**: the pipeline *Measured voxels → Mask → Symmetry average → Views → NEBULA3D*, followed by the Symmetry and Mask controls. Active stages are highlighted, and clicking a stage jumps to its controls.
  - **Export**: the last step: I(Q), the volume reduced to 1-D, and handing the processed volume to NEBULA3D.
- **Workspace**: a header row with the click mode, the display controls (see [Display](#display)) and the layout, then four views: HK, HL, KL, and the 3-D view, I(Q) or the line cut, which share the fourth place (the *3D | I(Q) | Cut* switch in its header picks one). On narrower windows the header row drops its button labels and colorbar (hover a button for its name), and wraps onto a second row when it still does not fit.
- **Panel footer**: copyright, the license and a link to this documentation. The ⓘ popover links to each guide.

The layout follows the screen:

- **Narrow screens** (below 1000 px, such as phones and tablets held upright): the views come first, one under the other (two per row from 700 px), each at most as tall as the screen, and the panel follows them; the page scrolls. The layout buttons are hidden, since every view is shown.
- **Short windows**: the panel stays beside the views, which keep a usable height; the views scroll.
- **Touch screens**: buttons, fields and sliders are 36–44 px tall, and two fingers pinch and pan a slice (see [Clicking a slice](#clicking-a-slice)).
- **Large screens at 100% scaling** (more than 2560×1440, such as 4K and 5K monitors): the whole interface, plot text included, is drawn 1.25–2× larger, so it is not tiny. Browser zoom (⌘ + / Ctrl +) adjusts it further.

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
| **Cut** | — | draws a [line cut](#line-cuts); a drag from an end of the cut moves that end |

In every mode, double-click (or *Reset zoom*) returns to the full view. On a touch screen, two fingers zoom by pinching and pan by moving together.

## Line cuts

A line cut is a 1-D profile through the volume along any line in 3-D: drawn in one of the slices (*3-D volume → 2-D slice → 1-D cut*), or typed, like one along (H, H, H). Each point averages a rod of voxels around the line, W across in every direction, not the slice's slab. A cut drawn in a slice lies in its plane, so it follows the slice as you move it; the slice's thickness does not change the cut.

- **Drawing**: choose **Cut** in the workspace header (or *Draw a cut* in the empty cut view) and drag across a slice. The ends snap to voxel centres; hold Shift to keep the cut along a lattice direction, such as (H, 0), (H, H) or (H, 2H). Drag an end to change the cut, or drag elsewhere to draw a new one. The first cut shows in the fourth place; after that the fourth place stays as you set it, so you can change a cut while watching the 3-D view. The cut shows on a slice it lies in as a line between two handles inside a dashed band, where its rod meets the slice's plane; on the other slices, where its rod passes through them: a dashed ellipse, its section, around a dot on its line; and in the 3-D view as a line through the slices inside the translucent rod of voxels it averages.
- **Typed ends** (footer): *From* and *To* take three coordinates. When they share one, like `-3, 0, 1` and `3, 0, 1`, the cut lies in the slice through it, which moves there: for these, the HK slice moves to L = 1 (the cut's current slice is kept when the ends fit it). When they share none, like `-2, -2, -2` and `2, 2, 2`, the cut is free in 3-D: here along (H, H, H), and the slices stay where they are.
- **Width W**: the diameter of the rod around the line, the same in every direction across it: in the slice's plane and out of it. It is in Å⁻¹ with a lattice (with 2π, as for I(Q)) or in the axes' unit otherwise. The slider covers half a voxel to 30 voxels on a log scale; empty means three voxels.
- **Points**: along the axis the cut changes most along (H for a cut along (H, 0.5H+1, 0) or (H, H, H)), every bin width of that axis by default; *Step* in the options changes it. Each point averages the voxels whose centres lie within W/2 of the line, in 3-D, and within half a step of the point, pooled with their symmetry equivalents as in the slices. See [Method → Line cuts](METHOD.md#line-cuts).
- **Plot**: as for I(Q): *Lin | Log*, hover to read the point's coordinates and I ± σ with the voxels pooled for each dataset (the point is marked on the slice), drag a box to zoom (a flat drag zooms along the cut only), and double-click for the full range. The options also hide the error band. σ comes from the file's uncertainties, read the first time a cut or I(Q) needs them.
- **Comparing**: both datasets are cut along the same line with their own grid, symmetry and mask, and drawn as two curves; *A / Split / B* chooses which are shown.
- **Download** (the arrow in the header) saves the points in a file named after the path, like `<file>_cut_H_0_1.csv` (`<A>_vs_<B>_cut_…` when comparing), as
  - **CSV**: a header row (`H,K,L,I,sigma,voxels`, or `I_A`, `sigma_A`, … for each dataset when comparing), then one row per point with its three coordinates and I, σ and voxels, and empty cells where there is no value; for spreadsheets and `pandas.read_csv`.
  - **Text (.dat)**: a commented header with the ends, slice, rod diameter, step, symmetry and mask, then one row per point with its position along the axis, its three coordinates, and I, σ and voxels for each dataset, `nan` where there is no value; for `numpy.loadtxt`.

  *Save PNG* exports the plot with a title and a legend.

## 3-D view

A transparent isosurface of the binned, symmetrized and masked volume, with the current slices as planes, clipped to the view range. The footer sets the isosurface level (log slider or typed; empty returns to the automatic level) and the surface and slice opacity. The options button sets the grid (about 64, 100 or 150 blocks per axis) and hides the slices. The [line cut](#line-cuts), when there is one, shows as a line through the slices inside a translucent tube, the rod of voxels it averages. The line is solid where it is in view and faint where a slice hides it, with a white collar where it crosses another slice. Drag to rotate, scroll to zoom, right-drag to pan, and use ↺ to reset the camera.

## I(Q)

I(Q) reduces the masked, symmetrized volume to one dimension: the mean intensity in each shell of |Q|, over the part of the shell that has data. Open it with *I(Q)* in the header of the 3-D view (the two share the fourth place), *Show I(Q)* in the **Export** section, or *Compute I(Q)* in its empty view. It is computed only once asked for, not when a file opens; after that, it is recomputed while shown when the symmetry, the mask or its settings change.

- **Normalization**: unmeasured and masked voxels are left out, not counted as zero, so gaps in coverage do not lower I(Q). Symmetry-equivalent voxels are pooled as in the slices, and each orbit counts with its multiplicity, so a shell is not biased toward the directions that were measured best. See [Method → I(Q)](METHOD.md#powder-average-iq).
- **Q bins** (footer): the ΔQ field takes Mantid `Rebin` parameters. `0.05` gives shells of 0.05 Å⁻¹; a negative step gives logarithmic shells, so `-0.01` means ΔQ/Q = 1% (suited to the constant Δd/d resolution of time-of-flight instruments); `0.5, 0.02, 3, 0.05, 10` gives 0.02 Å⁻¹ shells from 0.5 to 3 Å⁻¹ and 0.05 Å⁻¹ shells from 3 to 10 Å⁻¹, and any mix of ranges works. As in Mantid, a range ends at its boundary with a last bin of 0.25 to 1.25 steps. Empty means shells of the shortest bin step in |Q|.
- **ΔQ slider**: changes a single step on a log scale, from a tenth of the shortest bin step to 20 times it (or ΔQ/Q from 0.1% to 20% for logarithmic shells), and recomputes I(Q) as it moves. Typed values apply on Enter.
- **Q range**: Q min (empty for 0; the shortest bin step for logarithmic shells) and Q max (empty for all the data) bound a single step; ranges typed in the ΔQ field set their own, and the fields are then disabled. *Lin | Log* in the view's header switches the intensity scale. |Q| is in Å⁻¹ with 2π, from the file's UB matrix or cell; axes already in Å⁻¹ (Mantid's Q frames) need no cell.
- **Options** (the sliders button): *Split voxels* shares each voxel between the shells it overlaps by dividing it into 2³ or 3³ sub-cells, each binned by its own |Q|; *Centres* bins whole voxels, which is faster but aliases when ΔQ is close to the voxel size. The error band (±σ) and the shell coverage (the dashed line, right axis) can be hidden.
- **Uncertainties**: σ is propagated from the file's `errors_squared` (or `errors`), read the first time I(Q) is computed. Without them, the curve has no error band and the text file has `nan` in its σ column.
- **Coverage**: the fraction of each shell's volume that has data, after symmetry. Where it falls (beyond the measured region, or where shells leave the grid), I(Q) rests on few voxels.
- **Plot**: hover to read Q, d = 2π/Q, and I ± σ with the coverage of each dataset. Drag a box to zoom (a flat drag zooms Q only), click in *Zoom* mode to zoom 2×, drag in *Move* mode to pan, and double-click (or *Reset zoom*) for the full range.
- **Comparing**: both datasets are reduced on the same shells, each with its own grid, cell and mask, and drawn as two curves; *A / Split / B* chooses which are shown.
- **Download** (in the view's header or the Export section) saves a text file, `<file>_IQ.dat` (when comparing, the two names joined with their shared words once: `demo_300K` and `demo_10K` give `demo_300K_vs_10K_IQ.dat`): a commented header with the file, cell, shells (also as Rebin parameters, to reproduce them in Mantid), symmetry, mask and normalization, then one row per shell with Q (the shell centre, the midpoint of its edges), I, σ, coverage and the number of voxels with data, for A and B when comparing. `numpy.loadtxt` reads it directly. *Save PNG* exports the plot at 3× with a title and a legend.

## Comparing datasets

Up to four datasets can be compared, such as one crystal at four temperatures. *Compare…* in the top bar, or dropping a file on it, opens another dataset next to the first (A), as B, then C and D.

- **Two datasets** split every slice along its diagonal, from the top-left to the bottom-right corner of the view: A fills the lower-left half and B the upper-right half.
- **Three or four** share every slice in quadrants about the center of the view: A lower left, B lower right, C upper left and D upper right. With three, the fourth quadrant stays empty. With a symmetric pattern centered in the view, as for a Laue class with mirror planes, the quadrants show equivalent regions.

White gaps mark the cuts, the other datasets' parts are hatched where they have no data, and tags in the corners name the files. The quadrants follow the visible window, so zooming or panning moves the center they meet at.

Each dataset gets its own chip in the top bar, after A and the split icon, with its facts. While it loads, and while its mask is built, its chip shows the step, and the status in the top bar the progress; an amber *!* marks warnings (hover to read them).

- **Shared**: slice positions and thickness, zoom and pan, symmetry, mask settings and the color scale apply to every dataset, so their parts are directly comparable. *Auto range* pools them all.
- **A / Split / B / C / D** in the workspace header shows one dataset over the whole view, or the split.
- **Hover** reads every dataset at the cursor, the one under it first. View headers show the coverage of each.
- **I(Q) and line cuts** draw one curve per dataset (A blue, B amber, C green, D violet), and the saved files have columns for each.
- **Processing**: each dataset has its own worker and mask, built with the same parameters on its own data. Symmetry operations are mapped onto each grid; if they do not fit one, that dataset is used as measured and its chip says so.
- **Axes**: all datasets are drawn on A's axes and lattice geometry. Files should share the same axes (for example all `[H,0,0]`, `[0,K,0]`, `[0,0,L]`); the grids may differ. A dataset with A's axes in another order, such as a NEBULA3D volume (L, K, H) next to the Mantid file it came from (K, L, H), is read in A's order. A chip warns when its axis names differ from A's.
- **3-D view**: the isosurface is A's; the slice planes show the same split as the views.
- **Replace**: click a chip to swap in another file. **Remove** (× on the chip, which also cancels a loading dataset) frees its letter for the next *Compare…*; the others keep theirs. Opening a new file as A removes them all.

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
- **Remembered settings**: the colormap, scale, click mode, layout, the view in the fourth place, the I(Q) and line cut options, panel state and folded sections are remembered per browser.
