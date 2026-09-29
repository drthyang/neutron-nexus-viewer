# Method

[← README](../README.md) · [User guide](USER_GUIDE.md) · [Development](DEVELOPMENT.md)

The pipeline is *measured voxels → mask → symmetry average → slices, line cuts, 3-D view and I(Q)*. Masked and non-finite voxels are held as NaN. Voxels are stored as float32, and sums are accumulated in float64. There is no smoothing or interpolation at any stage.

## Slices

Each output pixel is the **equal-weight mean of the unique finite, unmasked source voxels** in the slab and all their symmetry equivalents. Symmetry orbits are either identical or disjoint, so the pooled set is the union of the distinct orbits met in the slab column. Each voxel is therefore counted once, including voxels that lie on symmetry elements (their images repeat |stabilizer| times, and sums are divided by that number). Pixels with no contribution are grey.

A slab is the set of bins whose centers lie within thickness/2 (+10⁻⁵) of the slab center.

## Line cuts

A line cut runs from a to b in the plane of a slice, within that slice's slab. It is sampled along its dominant axis, the in-plane axis it changes more along: point k is centred at a + k·step on that axis (the step defaults to the axis's bin width), and the last point is the one nearest b. A voxel of the slab belongs to point k when its centre lies within W/2 of the line and projects onto the line within step/2 of the point. Distances are measured in the plane's Cartesian metric, the same one the slices are drawn with (lengths |b\*ᵢ| per r.l.u. and the angle between the axes from G\*, with the nominal angles when those are drawn), so the band is the rectangle shown on the slice. W is given in Å⁻¹ with 2π on a lattice plane, and in the axes' unit otherwise.

Each point is then pooled exactly like a slice pixel: the **equal-weight mean of the unique finite, unmasked voxels in the union of the symmetry orbits** of its voxels, each orbit identified by its smallest in-range flat index and each distinct voxel counted once. The pooled voxels are independent measurements, so σ = √(Σσᵥ²)/N over the N voxels pooled, with σᵥ² from `errors_squared` (or `errors`). Drawn ends snap to voxel centres, so a cut along an axis between two voxel centres, with the default step, samples each voxel column once. A voxel whose centre projects exactly halfway between two points (as happens on diagonal cuts) goes to the upper one. With a second dataset, B is cut along the same line in A's display coordinates, on its own grid.

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

The volume is binned by an odd factor chosen so no axis has more than about 64, 100 or 150 blocks. Blocks are aligned to the origin: a block is centered on an origin bin, or has an edge on an origin edge. The same operations therefore map blocks onto blocks. Block sums and counts are pooled over each orbit, and the isosurface of the block means is drawn with surface nets. The default level is the 99.5th percentile of the positive block means. The line cut is drawn as a rod along its line, faint where a slice hides it, and the box it averages: its band across the line, through the slab of its slice.

## Powder average I(Q)

I(Q) is the spherical average of the masked, symmetrized volume. A voxel holds the intensity per unit reciprocal volume in its cell, so the average over the shell Q ≤ |Q| < Q + ΔQ is the **volume-weighted mean of the voxels in the shell that have data**:

I(Q) = Σ_v f_v I_v / Σ_v f_v, with f_v the part of voxel v inside the shell.

- **Coverage.** Equivalently, I(Q) is the intensity integrated over the part of the shell with data, divided by the volume of that part, not by the volume of the whole shell (their ratio is the reported coverage, so I·coverage is the integral over the shell divided by the whole shell's volume). Unmeasured and masked voxels are left out of both sums, not counted as zero, so a partly measured shell has the same expected intensity as a full one, only a larger uncertainty. Summing the shell instead would make I(Q) grow with Q² and with the coverage; dividing by the full shell volume would make it follow the coverage.
- **Symmetry.** An orbit's intensity is the equal-weight mean of its measured members, as in the slices, and every distinct member of the orbit counts with its own |Q| and volume, measured or not, inside the grid or beyond it. Orbits are therefore weighted by their multiplicity. Pooling only the measured voxels would weight each orbit by how many of its members happened to be measured, biasing the average toward the directions measured best: in a shell of three orbits of two voxels under mmm, where one orbit (intensity 100) is measured once and the others (intensity 1) fully, the shell average is (2·100 + 4·1)/6 = 34, where pooling the five measured voxels gives 20.8. Counting members beyond the grid also makes I(Q) independent of the shape of the grid's box.
- **Equal volumes, not inverse variances.** Weights are voxel volumes, as a spherical average requires. Inverse-variance weights would favour low-intensity voxels (σ² grows with the counts) and underestimate shells that contain Bragg peaks.
- **Shared voxels.** Each voxel is split into n³ sub-cells (2³ by default), each binned by its own |Q|, so a voxel is shared between the shells it overlaps. Binning whole voxels by their centres aliases when ΔQ is close to the voxel size: on a full 41³ grid with ΔQ equal to the shortest bin step, the coverage of the shells below 1.8 Å⁻¹ ranges from 0.64 to 1.31 with centres, and stays within 3% of 1 with 2³ sub-cells (the shell at the origin aside).
- **|Q|** = 2π √(hᵀG\*h), with h = Σ_d x_d b_d from the display coordinates x and the axes' basis vectors b (so `[H,H,0]` axes work), and G\* from the file's UB matrix or cell, not the nominal angles used for drawing. Axes in Å⁻¹ (Mantid's Q frames) are taken as Cartesian, with 2π already included.
- **Uncertainties** come from `errors_squared` (Mantid) or `errors` (NeXus standard deviations, squared on reading), read from the file when I(Q) is first computed. The members of an orbit share one mean and are fully correlated, so σ²(Q) = Σ_o c_o² σ_o² / (Σ_o c_o)², with σ_o² = Σ_m σ_m² / n_o² over the orbit's n_o measured members and c_o the orbit's weight in the shell. A shell with a voxel of unknown variance has σ = NaN.
- **Coverage** is reported per shell as the fraction of its volume with data after symmetry, Σ_v f_v V_voxel / (4π/3 ((Q+ΔQ)³ − Q³)). It is close to 1 for a measured shell inside the grid, and falls where shells leave the grid or the measured region, where I(Q) rests on few voxels.

**Shells** are given as Mantid `Rebin` parameters: one step ΔQ between Q min and Q max, or ranges Q₁, Δ₁, Q₂, Δ₂, Q₃ … each with its own step, where a negative step −r makes logarithmic shells, each edge (1 + r) times the previous (ΔQ/Q = r). The edges follow Mantid's rule: a range steps while the next edge plus a quarter step stays within its boundary, then ends at the boundary, so its last shell is 0.25 to 1.25 steps wide. Uniform edges are computed as Q₁ + kΔ, without accumulating rounding. By default ΔQ is the shortest bin step in |Q| (to one significant figure) and the shells run from 0 to the farthest grid corner; logarithmic shells start by default at the shortest bin step. Trailing shells without data are dropped. The coverage of each shell uses its own volume, 4π/3 (Q_hi³ − Q_lo³), and Q is reported at the shell centre, the midpoint of its edges. I(Q) is an angular average; a powder pattern per unit |Q| corresponds to 4πQ²·I(Q). With a second dataset, B is reduced on the same shells with its own grid, cell and mask.

## Comparing two datasets

A second dataset (B) is read by its own worker and processed independently with the same settings: the symmetry operations are converted to index maps on B's grid, and the mask is built from B's own data with the same erosion radius and outlier cut. Nothing is interpolated between the two grids. Each view draws both slices on A's axes and lattice geometry, clipped to the two triangles on either side of the view's diagonal, from (u₀, v₁) to (u₁, v₀). A point is in B's half when (u − u₀)/(u₁ − u₀) + (v − v₀)/(v₁ − v₀) > 1.

## Export for NEBULA3D

The export writes the input NEBULA3D's 3D-ΔPDF pipeline expects: `/entry/data` and `/entry/mask` (1 = valid) in (H, K, L) C order, bin-centre axes, and `ub_matrix` = 2π × the file's orientation matrix, or 2π × the Cholesky factor of G\* when the file has only a cell (only the metric matters for \|Q\| and the real-space axes). Each axis is padded to a grid symmetric about 0: odd with a bin centred at 0, or even when a bin edge lies at 0. Every voxel of the padded grid receives the equal-weight mean of the valid source voxels in its orbit, so the exported volume agrees with one-bin slices voxel by voxel (a test checks this for several Laue classes) and symmetry fills the padding. Orbits without a valid voxel are written as 0 with mask 0, NEBULA3D's convention for holes it backfills. Loaded in NEBULA3D, the example gives \|Q(100)\| = 1.7274 Å⁻¹ and \|Q(001)\| = 0.9240 Å⁻¹, as the cell requires, with equal intensity at 6/mmm-equivalent peaks.

## Geometry

The reciprocal metric is G\* = (UB)ᵀ·UB, with no 2π, taken from the file's UB matrix, or computed from `unit_cell_*` when there is no UB. The length of each axis per r.l.u. and the angle between two axes follow from G\* and the axes' HKL basis vectors (parsed from names such as `[H,H,0]`). Drawing is an affine map of the pixel grid, so bins keep their exact shape.

## Validation

- **Reference implementation**: with no symmetry or with −1, the slice engine reproduces `average_slab()` from the original Python viewer. The test suite checks it against fixtures generated by that function. On a real 401³ volume, counts and empty pixels matched exactly, and values agreed to within float32 precision (at most 3×10⁻⁵ relative).
- **Every Laue class** is checked against a direct set-based implementation on random volumes with gaps.
- **Line cuts** are checked against a direct set-based implementation for every Laue class, on axis-aligned, diagonal and oblique (hexagonal) lines, with σ and masks; an axis-aligned cut one voxel wide reproduces its row of voxels.
- **Masks**: edge erosion is checked against a brute-force box search, and the outlier cut against an injected spike.
- **Surface nets** is checked to produce a closed, consistently oriented sphere.
- **I(Q)** is checked against a direct set-based implementation (orbits as sets of distinct images, in the grid or not) for every Laue class, with and without split voxels, a mask and variances, and for logarithmic and piecewise shells; the shell edges are checked against Mantid's Rebin rules; against the three-orbit shell above; and on a smooth isotropic intensity, which it recovers within 2% at every shell inside the grid, also with 46% of the voxels removed.

Timing for a 401³ volume in a browser: loading about 1.5 s; a 3-bin slice with 6/mmm about 40 ms; a 41-bin slice about 0.5 s; building the mask about 5 s. I(Q) with 6/mmm and 2³ sub-cells takes about 8 s with uniform shells and 7 s with 1% logarithmic ones (measured in Node, with 70% of the voxels measured); the example takes under 0.1 s.
