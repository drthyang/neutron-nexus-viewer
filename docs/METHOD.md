# Method

[← README](../README.md) · [User guide](USER_GUIDE.md) · [Development](DEVELOPMENT.md)

The pipeline is *measured voxels → mask → symmetry average → slices and 3-D view*. Masked and non-finite voxels are held as NaN. Voxels are stored as float32, and sums are accumulated in float64. There is no smoothing or interpolation at any stage.

## Slices

Each output pixel is the **equal-weight mean of the unique finite, unmasked source voxels** in the slab and all their symmetry equivalents. Symmetry orbits are either identical or disjoint, so the pooled set is the union of the distinct orbits met in the slab column. Each voxel is therefore counted once, including voxels that lie on symmetry elements (their images repeat |stabilizer| times, and sums are divided by that number). Pixels with no contribution are grey.

A slab is the set of bins whose centers lie within thickness/2 (+10⁻⁵) of the slab center.

## Symmetry averaging

Operations are integer 3×3 matrices R acting on reflections, h′ = R·h. Real-space operations W are converted as R = Wᵀ. The generators are closed into a group (at most 48 operations), and two checks guard against mistakes:

- **Grid**: every operation, expressed in the HKL basis of the display axes, must send bin centers onto bin centers. It then becomes an exact integer map i′ = M·i + t on bin indices. Coupled axes need equal bin widths, and rotations other than inversion need a bin centered at the origin.
- **Metric**: the largest change of the reciprocal metric under the operations, max |RᵀG\*R − G\*| / max |G\*|, is reported. Values above 2% trigger a warning. A pseudo-symmetric cell shows a small but non-zero change; for example, 6/mmm changes a slightly distorted hexagonal cell by 0.8%.

Laue-class presets are given as generators in hexagonal axes for the trigonal and hexagonal classes: 1, −1, 2/m (b or c unique), mmm, 4/m, 4/mmm, −3, −3m1, −31m, 6/m, 6/mmm, m−3 and m−3m.

## Masking detector-edge artifacts

Detector edges, and the weak normalization there, leave spuriously high values along the boundaries of the measured region in reciprocal space. The reduced file no longer records which detector pixel a voxel came from, so both masks act on the unsymmetrized volume:

- **Edge erosion** removes measured voxels within *r* voxels (box distance) of an unmeasured voxel. It is a separable dilation of the unmeasured set along each axis.
- **Outlier cut** visits each symmetry orbit with at least 3 valid voxels. Voxels above median + *k*·1.4826·MAD of their orbit are removed. The outlier decisions are made on the voxels left after erosion.

In a 401³ Fe₃Ge₂ volume (90 K), the 99th percentile of voxels within one voxel of an edge was 4266, against 36.5 for interior voxels. Erosion by 1 voxel removed 15.8% of measured voxels, and erosion by 2 removed 30.6%. With *r* = 2 and a 5σ outlier cut under 6/mmm, the HK plane at L = 0 had 120 pixels above 1000 instead of 1476. Its coverage fell only from 73.9% to 71.7%, because symmetry fills most gaps.

The cleaner fix is upstream: mask detector-edge pixels (for example with Mantid's `MaskBTP`) before converting to MD and normalizing.

## 3-D view

The volume is binned by an odd factor chosen so no axis has more than about 64, 100 or 150 blocks. Blocks are aligned to the origin: a block is centered on an origin bin, or has an edge on an origin edge. The same operations therefore map blocks onto blocks. Block sums and counts are pooled over each orbit, and the isosurface of the block means is drawn with surface nets. The default level is the 99.5th percentile of the positive block means.

## Comparing two datasets

A second dataset (B) is read by its own worker and processed independently with the same settings: the symmetry operations are converted to index maps on B's grid, and the mask is built from B's own data with the same erosion radius and outlier cut. Nothing is interpolated between the two grids. Each view draws both slices on A's axes and lattice geometry, clipped to the two triangles on either side of the view's diagonal, from (u₀, v₁) to (u₁, v₀). A point is in B's half when (u − u₀)/(u₁ − u₀) + (v − v₀)/(v₁ − v₀) > 1.

## Export for NEBULA3D

The export writes the input NEBULA3D's 3D-ΔPDF pipeline expects: `/entry/data` and `/entry/mask` (1 = valid) in (H, K, L) C order, bin-centre axes, and `ub_matrix` = 2π × the file's orientation matrix, or 2π × the Cholesky factor of G\* when the file has only a cell (only the metric matters for \|Q\| and the real-space axes). Each axis is padded to a grid symmetric about 0: odd with a bin centred at 0, or even when a bin edge lies at 0. Every voxel of the padded grid receives the equal-weight mean of the valid source voxels in its orbit, so the exported volume agrees with one-bin slices voxel by voxel (a test checks this for several Laue classes) and symmetry fills the padding. Orbits without a valid voxel are written as 0 with mask 0, NEBULA3D's convention for holes it backfills. Loaded in NEBULA3D, the example gives \|Q(100)\| = 1.7274 Å⁻¹ and \|Q(001)\| = 0.9240 Å⁻¹, as the cell requires, with equal intensity at 6/mmm-equivalent peaks.

## Geometry

The reciprocal metric is G\* = (UB)ᵀ·UB, with no 2π, taken from the file's UB matrix, or computed from `unit_cell_*` when there is no UB. The length of each axis per r.l.u. and the angle between two axes follow from G\* and the axes' HKL basis vectors (parsed from names such as `[H,H,0]`). Drawing is an affine map of the pixel grid, so bins keep their exact shape.

## Validation

- **Reference implementation**: with no symmetry or with −1, the slice engine reproduces `average_slab()` from the original Python viewer. The test suite checks it against fixtures generated by that function. On a real 401³ volume, counts and empty pixels matched exactly, and values agreed to within float32 precision (at most 3×10⁻⁵ relative).
- **Every Laue class** is checked against a direct set-based implementation on random volumes with gaps.
- **Masks**: edge erosion is checked against a brute-force box search, and the outlier cut against an injected spike.
- **Surface nets** is checked to produce a closed, consistently oriented sphere.

Timing for a 401³ volume in a browser: loading about 1.5 s; a 3-bin slice with 6/mmm about 40 ms; a 41-bin slice about 0.5 s; building the mask about 5 s.
