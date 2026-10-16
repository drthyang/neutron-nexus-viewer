# Features

[← README](../README.md) · [User guide](USER_GUIDE.md) · [Method](METHOD.md) · [Agent tools](AGENT.md) · [Development](DEVELOPMENT.md)

What the viewer is for, everything it does with links to the details, and a first session with it.

## Goals

- **The data as measured, in its true geometry.** Slices are drawn with the
  reciprocal metric from the file's UB matrix, so oblique axes meet at their real
  angle and equal lengths in r.l.u. look equal. Slices and cuts are never smoothed or
  interpolated; only reindexing resamples the volume, which smooths it (peaks about 9% lower
  on the example).
- **Processing you can check.** Symmetry averaging maps bins onto bins exactly and
  warns when the operations do not fit the cell. Masks can be previewed voxel by
  voxel and exported. The slice engine is tested against the reference Python
  implementation. See [Method](METHOD.md).
- **Fast, private, nothing to install.** A static page reads the file locally with
  [h5wasm](https://github.com/usnistgov/h5wasm): a 401³ volume opens in about 2 s,
  and slices update in tens of milliseconds.

## Features

| Area | What you get | Read more |
| --- | --- | --- |
| **Slices** | HK, HL and KL cuts with live slab center and thickness. Click to move the other slices through a point, drag a box to zoom, or drag to pan (pinch on a touch screen); quad, focus and single layouts | [Slice views](USER_GUIDE.md#slice-views) |
| **Line cuts** | Drag a line across a slice, or type its ends, for a 1-D profile of the volume along any direction in 3-D, such as (H, H, H), averaged over a rod of adjustable diameter around the line, with symmetry pooling, σ and both datasets; saved as CSV or text | [Line cuts](USER_GUIDE.md#line-cuts) · [Method](METHOD.md#line-cuts) |
| **Symmetry averaging** | Any Laue class or your own operations, closed into a group and applied exactly on the bin grid, with a check that they fit the cell | [Symmetry](METHOD.md#symmetry-averaging) |
| **Artifact masking** | Removes detector-edge voxels and symmetry outliers before averaging; preview what is removed and export the mask as a NumPy array | [Masking](METHOD.md#masking-detector-edge-artifacts) |
| **UB refinement** | Refine the UB against the volume's own Bragg peaks, with the cell constrained to a crystal system (cubic to triclinic), standard uncertainties and the unconstrained fit to compare; save it as an ISAW UB file or reindex with it | [Refining the UB](USER_GUIDE.md#refining-the-ub) · [Method](METHOD.md#refining-the-ub) |
| **Reindexing with a new UB** | Load an ISAW UB file (Mantid's `SaveIsawUB`), such as a refined UB or another cell setting, and resample the volume onto its HKL grid with σ propagated; opens for comparison and saves as a Mantid file | [Reindexing](USER_GUIDE.md#reindexing-with-a-new-ub) · [Method](METHOD.md#reindexing-with-a-new-ub) |
| **Comparing up to four datasets** | Other files, such as other temperatures: a second one splits every slice along its diagonal, three or four share it in quadrants, with shared positions, processing and color scale, or a color range per dataset for data on different scales, such as X-ray next to neutron; hover reads every value, and I(Q) and line cuts overlay them | [Comparing](USER_GUIDE.md#comparing-datasets) |
| **3-D view** | A transparent isosurface of the processed volume, with the current slices as planes and the line cut, inside the rod it averages, passing through them | [3-D view](USER_GUIDE.md#3-d-view) |
| **I(Q)** | The processed volume reduced to 1-D, on request: the mean intensity per \|Q\| shell over the voxels with data, with symmetry orbits weighted by multiplicity, split voxels, propagated σ and shell coverage; saved as text | [I(Q)](USER_GUIDE.md#iq) · [Method](METHOD.md#powder-average-iq) |
| **Hand-off to NEBULA3D** | One click sends the masked, symmetrized volume to [NEBULA3D](https://github.com/drthyang/nebula3d)'s 3D-ΔPDF pipeline in a new tab, or saves it as a file | [Export](USER_GUIDE.md#export-for-nebula3d) |
| **Any screen** | Laptops to 4K and 5K monitors (drawn larger at 100% scaling), tablets and phones: views first on narrow screens, touch-sized controls and pinch zoom | [The screen](USER_GUIDE.md#the-screen) |
| **Rigaku XRD reduction** | Raw CrysAlisPro frames (`*.rod_img`) to an HKL volume in the browser: detector mask, Bragg-peak search, geometry refinement per run, solid-angle and polarization normalization, in any cell setting (a transformation matrix) or on an open dataset's grid, optional background subtraction (measured, static or rotation-resolved, or an isotropic one such as a Kapton mount's estimated from the frames), using all CPU cores with the same result as one; opens for comparison and saves as a Mantid file that `LoadMD` reads | [Reducing frames](USER_GUIDE.md#reducing-rigaku-xrd-frames) · [Method](METHOD.md#rigaku-reduction) |
| **Files and links** | Mantid `MDHistoWorkspace` (`SaveMD`), any 3-D `NXdata`, and NEBULA3D's processed volumes and 3D-ΔPDFs (drawn at the real-space cell angles); local files or links that open a file, a comparison, symmetry and mask; PNG export with a colorbar | [Supported files](USER_GUIDE.md#supported-files) |
| **Agent tools** | The analysis as tools for AI agents over MCP (Claude Code, Claude Desktop, any MCP client): open a volume, set symmetry and a mask, look at slices, take line cuts and I(Q), read the intensity at NEXPLAN's Bragg nodes, export for NEBULA3D, with the viewer's own code | [Agent tools](AGENT.md) |

## A first session

1. Open the [live app](https://drthyang.github.io/neutron-nexus-viewer/) and choose a `.nxs` file, or drop one on the page.
2. Pick a Laue class under **Processing → Symmetry averaging**, and apply a
   **Mask** if detector edges show up as bright rims.
3. Click a slice to move the other two through that point. Switch the header to
   **Zoom** or **Move** to explore, and **Focus** to show one view large.
4. For a 1-D profile, switch the header to **Cut** and drag across a slice: the cut
   appears in the fourth view, and follows the slice as you move it.
5. To compare, click **Compare…** next to the file in the top bar and open a second file: each slice is split along its diagonal. Add a third and fourth file to share each slice in quadrants.
6. For a 3D-ΔPDF, choose **Open in NEBULA3D** in the **Export** section:
   [NEBULA3D](https://drthyang.github.io/nebula3d/) opens with the volume loaded.
