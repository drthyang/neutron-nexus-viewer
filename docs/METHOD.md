# Method

[← README](../README.md) · [User guide](USER_GUIDE.md) · [Development](DEVELOPMENT.md)

The pipeline is *measured voxels → mask → symmetry average → slices, line cuts, 3-D view and I(Q)*. Masked and non-finite voxels are held as NaN. Voxels are stored as float32, and sums are accumulated in float64. There is no smoothing or interpolation at any stage.

## Slices

Each output pixel is the **equal-weight mean of the unique finite, unmasked source voxels** in the slab and all their symmetry equivalents. Symmetry orbits are either identical or disjoint, so the pooled set is the union of the distinct orbits met in the slab column. Each voxel is therefore counted once, including voxels that lie on symmetry elements (their images repeat |stabilizer| times, and sums are divided by that number). Pixels with no contribution are grey.

A slab is the set of bins whose centers lie within thickness/2 (+10⁻⁵) of the slab center.

## Line cuts

A line cut runs from a to b, any line in 3-D: drawn in a slice, and so in that slice's plane, or typed. It is sampled along its dominant axis, the axis it changes most along: point k is centred at a + k·step on that axis (the step defaults to the axis's bin width), and the last point is the one nearest b. A voxel belongs to point k when its centre lies inside a **rod of diameter W** around the line, within W/2 of it in any direction across it, and projects onto the line within step/2 of the point. Distances are measured in reciprocal space, with the Cartesian basis the 3-D view uses (the reciprocal lattice vectors from G\*, with the nominal angles when those are drawn), so the rod is round whatever the lattice. On a slice whose plane holds the line, the band drawn is where the rod meets that plane; on a slice the line crosses at an angle θ to the plane's normal, the rod's section is an ellipse, W across the line's shadow and W / cos θ along it, drawn around the crossing point. The rod replaces the slice's slab: its depth out of the plane is W, whatever the slice's thickness. W is given in Å⁻¹ with 2π with a lattice, and in the axes' unit otherwise.

Each point is then pooled exactly like a slice pixel: the **equal-weight mean of the unique finite, unmasked voxels in the union of the symmetry orbits** of its voxels, each orbit identified by its smallest in-range flat index and each distinct voxel counted once. The pooled voxels are independent measurements, so σ = √(Σσᵥ²)/N over the N voxels pooled, with σᵥ² from `errors_squared` (or `errors`). Drawn ends snap to voxel centres, so a cut along an axis between two voxel centres, with the default step, samples each voxel column once. A voxel whose centre projects exactly halfway between two points (as happens on diagonal cuts) goes to the upper one. With other datasets, each is cut along the same line in A's display coordinates, on its own grid.

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

The volume is binned by an odd factor chosen so no axis has more than about 64, 100 or 150 blocks. Blocks are aligned to the origin: a block is centered on an origin bin, or has an edge on an origin edge. The same operations therefore map blocks onto blocks. Block sums and counts are pooled over each orbit, and the isosurface of the block means is drawn with surface nets. The default level is the 99.5th percentile of the positive block means. The line cut is drawn as a line, faint where a slice hides it, inside a translucent tube: the rod of voxels it averages.

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

**Shells** are given as Mantid `Rebin` parameters: one step ΔQ between Q min and Q max, or ranges Q₁, Δ₁, Q₂, Δ₂, Q₃ … each with its own step, where a negative step −r makes logarithmic shells, each edge (1 + r) times the previous (ΔQ/Q = r). The edges follow Mantid's rule: a range steps while the next edge plus a quarter step stays within its boundary, then ends at the boundary, so its last shell is 0.25 to 1.25 steps wide. Uniform edges are computed as Q₁ + kΔ, without accumulating rounding. By default ΔQ is the shortest bin step in |Q| (to one significant figure) and the shells run from 0 to the farthest grid corner; logarithmic shells start by default at the shortest bin step. Trailing shells without data are dropped. The coverage of each shell uses its own volume, 4π/3 (Q_hi³ − Q_lo³), and Q is reported at the shell centre, the midpoint of its edges. I(Q) is an angular average; a powder pattern per unit |Q| corresponds to 4πQ²·I(Q). With other datasets, each is reduced on the same shells with its own grid, cell and mask.

## Comparing datasets

Each further dataset (B, C, D) is read by its own worker and processed independently with the same settings. When its axes are A's in another order (matched by name), its volume is transposed to A's order as it is read, so all share the display axes: the symmetry operations are converted to index maps on its grid, and the mask is built from its own data with the same erosion radius and outlier cut. Nothing is interpolated between grids. Each view draws the slices on A's axes and lattice geometry, each clipped to its part of the visible window [u₀, u₁] × [v₀, v₁]:

- **Two datasets**: the two triangles on either side of the view's diagonal, from (u₀, v₁) to (u₁, v₀). A point is in B's half when (u − u₀)/(u₁ − u₀) + (v − v₀)/(v₁ − v₀) > 1.
- **Three or four**: the quadrants about the window's center (u_m, v_m) = ((u₀ + u₁)/2, (v₀ + v₁)/2), cut along the plot axes, so on oblique axes they are parallelograms. A point belongs to A when u < u_m and v < v_m, to B when u ≥ u_m and v < v_m, to C when u < u_m and v ≥ v_m, and to D otherwise.

All datasets are colored on one scale by default, and *Auto* takes its percentiles over the values of all their slices pooled. With own color ranges, dataset k is colored with its own limits and softening, and *Auto* takes them from its slices alone, as for a single dataset; the colormap and the scale function stay shared. Equal colors then mean equal intensity only within a dataset.

I(Q) and line cuts of every dataset use A's shells or A's line in display coordinates, each on its own grid.

## Export for NEBULA3D

The export writes the input NEBULA3D's 3D-ΔPDF pipeline expects: `/entry/data` and `/entry/mask` (1 = valid) in (H, K, L) C order, bin-centre axes, and `ub_matrix` = 2π × the file's orientation matrix, or 2π × the Cholesky factor of G\* when the file has only a cell (only the metric matters for \|Q\| and the real-space axes). Each axis is padded to a grid symmetric about 0: odd with a bin centred at 0, or even when a bin edge lies at 0. Every voxel of the padded grid receives the equal-weight mean of the valid source voxels in its orbit, so the exported volume agrees with one-bin slices voxel by voxel (a test checks this for several Laue classes) and symmetry fills the padding. Orbits without a valid voxel are written as 0 with mask 0, NEBULA3D's convention for holes it backfills. Loaded in NEBULA3D, the example gives \|Q(100)\| = 1.7274 Å⁻¹ and \|Q(001)\| = 0.9240 Å⁻¹, as the cell requires, with equal intensity at 6/mmm-equivalent peaks.

## Rigaku reduction

**Reduce Rigaku XRD…** turns the raw frames of a Rigaku Oxford Diffraction (CrysAlisPro) single-crystal experiment into an HKL volume, in a worker (`js/rigaku-worker.js`, `js/rigaku-format.js`, `js/rigaku-geometry.js`, `js/rigaku-reduce.js`). It reads every frame three times and never holds more than a few detector-sized arrays and the output grid.

**Inputs.** Files are found by name from the frames' experiment stem:
- the main-series frames `<stem>_<run>_<frame>.rod_img`; `pre_*` screening runs are listed but not used;
- `<stem>.par`, for the monochromator;
- `<stem>_cracker.par` and `expinfo/<stem>_crystal.ini`, for the orientation matrices and the Laue class;
- `expinfo/<stem>_datacoll.ini`, for the temperature.

**Frames.** The header ("OD SAPPHIRE", offsets as in FabIO and dxtbx) gives the scan, the goniometer angles and zero corrections, the exposure, the wavelength and the detector model. The TY6-compressed pixels are decoded as in dxtbx's `FormatROD`. Every decoded frame must reproduce the min, max, mean and standard deviation that the instrument software stored in its header, or the reduction stops. The stored counts carry no dark, flat-field or geometric correction (header `Unwarping` = 0, flood `NONE`).

**Geometry.** The conventions are those of dxtbx:
- *Laboratory frame:* the beam travels along −Z, +Y is along −ω, and rotations are right-handed.
- *Goniometer:* x_lab = R(ω) R(κ) R(φ) x_C, with ω and φ about (0, −1, 0) and κ about (0, −cos α, sin α).
- *Detector:* R_det = R(−Y, 2θ_arm) R(−X, d₂) R(Z, d₁); pixel (i, j) lies at R_det[(i − x₀)p, (j − y₀)p, −D].
- *Zero corrections:* the header's software zero corrections are added to the motor angles; on the XtaLAB mini II, ω gets +90°.
- *UB frame:* the CrysAlis UB includes the wavelength (|UB·h| = λ/d) and is expressed in its own frame e₁ = +Z, e₂ = +X, e₃ = +Y.
- *Indices:* hkl = UB⁻¹ M x_C, with M = [[0,0,1],[1,0,0],[0,1,0]] and x = s − s₀ (no 2π).

These choices were fixed empirically. For two XtaLAB mini II data sets (10 runs each, κ = 54°, ten φ settings), all 4,608 and 6,144 combinations of axis permutation, ω zero offset and rotation senses were scored by the fraction of independently found peaks that index to integers. Only this one indexes all runs: 99.9 % of 815 peaks on one data set, against a median of 1 % for the other combinations. It is exact up to Friedel inversion, which a centrosymmetric Laue class cannot distinguish.

**Detector mask.** Built from the per-pixel sum over all frames:
- the outer border;
- chip-boundary triplets, where a row or column deviates by more than 8 % from a 15-line median trend; on a 775 × 385 HyPix-3000, its known boundaries (columns 96, 193, 290, 387, 484, 581, 678 and row 192, each ± 1) are always included, because their contrast varies between data sets;
- the beamstop umbra: pixels below 20 % of a 41-pixel median level, in components of at least 50 px, dilated by 3;
- its penumbra: pixels connected to the umbra and below 75 % of a 61-pixel 90th-percentile level, at most 26 px from the umbra, dilated by 2;
- pixels that never counted.

**Bragg peaks.** A pixel is strong when it has ≥ 8 counts and ≥ b + 6√(b+1), where b is the pixel's mean over its run. Strong pixels are joined across their 4 neighbours and the same pixel in the next frame. Components with ≥ 4 voxels and ≥ 150 net counts are kept, at net-weighted centroids, with frame midpoints for the scan angle. Peaks touching the mask are left out.

**Refinement.** The starting UB is whichever of the CrysAlis matrices indexes most peaks (within 0.1). The cell is constrained by the crystal system of the Laue class, when the starting cell fits it. The fit is Levenberg–Marquardt with a soft-L1 loss on detector x, y (pixels) and scan angle (0.1° weighted as 1 pixel), in three stages:
- **L3:** orientation, cell, beam centre, in-plane detector rotation d₁, and the scan-axis zero.
- **L6:** adds the crystal's offset from the rotation centre and the κ zero.
- **L10:** the goniometer is fixed, and each run gets its own small orientation correction relative to the run with most peaks. A run whose angular rms stays above max(0.15°, 2 × median) gets a piecewise-linear drift with 6 knots instead.

d₂ stays at the header value: it is degenerate with the beam centre. The detector distance stays at the header's calibrated value too. The peaks fix only the ratio of cell to distance, so refining the distance moves the absolute cell scale without improving the fit. On one data set it changed the rms by 0.01 px and the HKL map by at most 0.006 r.l.u., but changed a by 0.12 %.

**Parallel workers.** With more than two CPU cores, the frames of each pass are spread over Web Workers (one per core, leaving one free). The result is bit-identical to a serial run:
- pass 1 sums integer counts, which is exact in any order;
- pass 2 gives each worker whole runs, so the peak search is unchanged;
- in pass 3 the workers compute each pixel's voxels and contributions. One thread adds them frame by frame, grouped by blocks of 4,096 voxels with a stable sort, so every voxel receives the same terms in the same order as in a serial run. Floating-point sums depend on that order, so the order is kept.

Where workers cannot start, the reduction runs serially with the same result.

**Gridding.** Each unmasked pixel of each frame is split into n equal sub-steps of the frame's rotation (default 5). The shutterless frame integrates continuously, so its counts are shared equally among the sub-steps, and each sub-sample goes to the voxel containing it. The crystal offset is evaluated at the frame midpoint. Per voxel:
- S = Σ f·c (counts);
- E2 = Σ f²·c, where f is a pixel-frame's share in the voxel, summed over its sub-samples before squaring, so split counts are not counted as independent;
- W = Σ f·t·w_p, the normalization weight;
- N, the number of pixel-frames.

The signal is S/W with errors² = E2/W². The weight w_p is either 1 (*exposure only*: counts per second per pixel), or ΔΩ_p·P_p/Ω_ref (*solid angle + polarization*), where:
- ΔΩ_p = p² cos α / r² is the pixel's solid angle;
- P_p is the polarization factor of the graphite-monochromated beam, [(1 − s_σ²) + cos²2θ_m (1 − s_π²)]/(1 + cos²2θ_m), with σ perpendicular to the monochromator plane given in the `.par`, or (1 + cos²2θ)/2 without a monochromator;
- Ω_ref = (p/D)².

Signal values are voxel averages of a continuous-scattering estimate, so diffuse scattering needs no Lorentz factor; Bragg-peak voxels are not integrated intensities. Background (air scatter, fluorescence), absorption and symmetry averaging are not applied. An absorption correction needs a crystal shape, which CrysAlis files usually lack.

**Grid and output.**
- *Indices:* the output cell is a transformation T of the refined cell. Each row of T gives an output basis vector in units of a, b, c (a′ᵢ = Σⱼ Tᵢⱼ aⱼ), and the output indices are T · (h, k, l). Entries may be fractions. The determinant must be positive, so the cell stays right-handed and its orientation is a proper rotation for Mantid. Presets:
  - identity;
  - 2 × 2 × 2, the doubled cell of many neutron reductions;
  - orthohexagonal (a′ = a, b′ = a + 2b, c′ = c), for hexagonal cells;
  - *Match dataset A*, n × identity with n from A's cell, and voxel centres on A's.

  The file stores the transformed cell and UB, and T in a log.
- *Extent:* the default range is everything the detector reaches.
- *File:* the result is written in the layout of Mantid's `SaveMD` (version 2): `signal`, `errors_squared`, `num_events` (pixel-frames) and `mask` in (L, K, H) order with `axes = D2:D1:D0`, HKL dimensions, and the oriented lattice and UB under `experiment0`. The reduction report is stored as a log.
- *Unmeasured voxels:* NaN signal and errors, with `num_events` = 0. Measured zeros stay 0.

**Validation**, on two XtaLAB mini II / HyPix-3000 data sets (7,342 frames):
- *Decoder:* every frame matches its header statistics, and frames are identical, pixel by pixel, to an independent Python decoder.
- *Geometry:* forward predictions agree with a Python implementation to 10⁻⁸ px.
- *Gridding:* with the same geometry and mask, the accumulators agree with the reference Python implementation voxel by voxel over 3.3 million voxels: S and W to 10⁻¹⁴, E2 to float32 precision, identical N and coverage.
- *Refinement:* the in-browser refinement reproduces the Python one at every level. Final rms is 0.361 / 0.394 px and 0.191° (Python 0.359 / 0.392 px, 0.188°), with cell a = 4.0082, c = 5.0156 Å (Python 4.0083 / 5.0156), both with the distance at the header value.
- *Volume:* Bragg integrated intensities agree within 1 % (median 0.998).
- *Mantid:* `LoadMD` (Mantid 6.16.1) reads the file with its dimensions, frame, lattice and logs. SliceViewer draws HK planes at 60°.
- *Parallel workers:* serial and parallel reductions give byte-identical accumulators and the same refined model. This was checked on simulated frames (tests, and in a browser with real Web Workers and nested workers) and on 736 and 3,608 real frames.
- *Speed:* 3,608 frames take about 40 s with 7 workers, against about 2 minutes serially, at about 1.5 GB.
- *Tests:* `tests/rigaku.test.js` checks each step on synthetic frames, including a simulated experiment whose refined geometry puts the Bragg peaks back on integer HKL. `tests/rigaku-local.test.js` checks every frame of a real experiment when `RIGAKU_DIR` is set.

## Geometry

The reciprocal metric is G\* = (UB)ᵀ·UB, with no 2π, taken from the file's UB matrix, or computed from `unit_cell_*` when there is no UB. The length of each axis per r.l.u. and the angle between two axes follow from G\* and the axes' HKL basis vectors (parsed from names such as `[H,H,0]`). Drawing is an affine map of the pixel grid, so bins keep their exact shape.

Real-space axes x, y, z in Å along a, b and c (a NEBULA3D 3D-ΔPDF, where the point (x, y, z) is at x·â + y·b̂ + z·ĉ) are drawn at the direct-cell angles: α between y and z, β between x and z, γ between x and y. The 3-D view places them along the unit cell vectors, and line-cut widths are in Å. A symmetry operation R given on h, k, l acts on fractional coordinates as W = Rᵀ, so on these axes as x′ = D Rᵀ D⁻¹ x with D = diag(\|a\|, \|b\|, \|c\|); without a cell, operations act on the axes directly.

## Validation

- **Reference implementation**: with no symmetry or with −1, the slice engine reproduces `average_slab()` from the original Python viewer. The test suite checks it against fixtures generated by that function. On a real 401³ volume, counts and empty pixels matched exactly, and values agreed to within float32 precision (at most 3×10⁻⁵ relative).
- **Every Laue class** is checked against a direct set-based implementation on random volumes with gaps.
- **Line cuts** are checked against a direct set-based implementation for every Laue class, on axis-aligned, diagonal and oblique (hexagonal) lines, in and out of a slice's plane, with σ and masks; a thin rod along an axis reproduces its row of voxels, and a wider one takes the same neighbours in and out of the plane.
- **Masks**: edge erosion is checked against a brute-force box search, and the outlier cut against an injected spike.
- **Surface nets** is checked to produce a closed, consistently oriented sphere.
- **I(Q)** is checked against a direct set-based implementation (orbits as sets of distinct images, in the grid or not) for every Laue class, with and without split voxels, a mask and variances, and for logarithmic and piecewise shells; the shell edges are checked against Mantid's Rebin rules; against the three-orbit shell above; and on a smooth isotropic intensity, which it recovers within 2% at every shell inside the grid, also with 46% of the voxels removed.

Timing for a 401³ volume in a browser: loading about 1.5 s; a 3-bin slice with 6/mmm about 40 ms; a 41-bin slice about 0.5 s; building the mask about 5 s. I(Q) with 6/mmm and 2³ sub-cells takes about 8 s with uniform shells and 7 s with 1% logarithmic ones (measured in Node, with 70% of the voxels measured); the example takes under 0.1 s.
