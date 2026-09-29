import assert from 'node:assert/strict';
import test from 'node:test';

import { parseBasis, reciprocalMetric } from '../js/nexus.js';
import { parseBins, powderAverage, powderPlan, qMetric, shellEdges } from '../js/powder.js';
import { IDENTITY_MAP } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps, PRESETS } from '../js/symmetry.js';

const group = (name) => closeGroup(parseOps(PRESETS.find(([n]) => n === name)[1]));
// An HKL axis with n bins of width w centred on the origin (or starting at `start`).
const axis = (letter, n, w, start = -n * w / 2) => ({
  label: letter, frame: 'HKL', units: 'r.l.u.', basis: parseBasis(`[${[...'HKL'].map((c) => (c === letter ? c : 0))}]`),
  edges: Float64Array.from({ length: n + 1 }, (_, i) => start + i * w),
});
const hexCell = { a: 4.2, b: 4.2, c: 6.8, alpha: 90, beta: 90, gamma: 120 };
const cubic = { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 };

function random(seed) {
  let s = seed;
  return () => ((s = (s * 16807) % 2147483647) / 2147483647);
}

function randomVolume(n, seed, gaps = 0.3) {
  const rand = random(seed);
  return Float32Array.from({ length: n }, () => (rand() < gaps ? NaN : rand() * 10));
}

/** A volume of display shape n with values f(x) at bin centres x (display coordinates). */
function fill(dims, f) {
  const n = dims.map((d) => d.edges.length - 1), volume = new Float32Array(n[0] * n[1] * n[2]);
  const c = (d, i) => (dims[d].edges[i] + dims[d].edges[i + 1]) / 2;
  for (let i2 = 0, p = 0; i2 < n[2]; i2++) for (let i1 = 0; i1 < n[1]; i1++) for (let i0 = 0; i0 < n[0]; i0++, p++) volume[p] = f([c(0, i0), c(1, i1), c(2, i2)]);
  return { volume, shape: [n[2], n[1], n[0]] };
}

/**
 * Direct implementation: orbits as sets of distinct images (in the grid or
 * not), each distinct image split into sub-cells binned by its own |Q|.
 */
function reference(volume, shape, plan, maps, mask, variance) {
  const n = [shape[2], shape[1], shape[0]], { G, c0, w, edges, count, split } = plan;
  const inGrid = (q) => q.every((x, i) => x >= 0 && x < n[i]);
  const flat = (q) => q[0] + n[0] * (q[1] + n[1] * q[2]);
  const qOf = (x) => Math.sqrt(Math.max(0, [0, 1, 2].reduce((s, i) => s + x[i] * [0, 1, 2].reduce((t, j) => t + G[3 * i + j] * x[j], 0), 0)));
  const weight = new Float64Array(count), sum = new Float64Array(count), spread = new Float64Array(count);
  const seen = new Set();
  for (let i2 = 0; i2 < n[2]; i2++) for (let i1 = 0; i1 < n[1]; i1++) for (let i0 = 0; i0 < n[0]; i0++) {
    const p = [i0, i1, i2];
    if (seen.has(flat(p))) continue;
    const members = new Map();
    for (const { M, t } of maps) {
      const q = [0, 1, 2].map((i) => M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + t[i]);
      members.set(q.join(), q);
    }
    const measured = [];
    for (const q of members.values()) {
      if (!inGrid(q)) continue;
      seen.add(flat(q));
      const v = volume[flat(q)];
      if (!Number.isNaN(v) && !(mask && mask[flat(q)])) measured.push(flat(q));
    }
    if (!measured.length) continue;
    const mean = measured.reduce((s, f) => s + volume[f], 0) / measured.length;
    const error = variance ? measured.reduce((s, f) => s + variance[f], 0) / measured.length ** 2 : 0;
    const shells = new Map();
    for (const q of members.values()) {
      for (let a = 0; a < split; a++) for (let b = 0; b < split; b++) for (let c = 0; c < split; c++) {
        const x = [a, b, c].map((k, i) => c0[i] + (q[i] + (k + 0.5) / split - 0.5) * w[i]);
        const Q = qOf(x), shell = edges.findIndex((e, s) => s < count && Q >= e && Q < edges[s + 1]);
        if (shell >= 0) shells.set(shell, (shells.get(shell) ?? 0) + 1 / split ** 3);
      }
    }
    for (const [shell, c] of shells) {
      weight[shell] += c;
      sum[shell] += c * mean;
      spread[shell] += c * c * error;
    }
  }
  return {
    intensity: Array.from(weight, (x, b) => (x ? sum[b] / x : NaN)),
    sigma: Array.from(weight, (x, b) => (x && variance ? Math.sqrt(spread[b]) / x : NaN)),
    voxels: weight,
  };
}

const close = (a, b, tol = 1e-9) => (Number.isNaN(a) ? Number.isNaN(b) : Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)));

test('every Laue class matches a direct set-based implementation', () => {
  // L has fewer bins than H and K, so the cubic classes send members out of the grid.
  const dims = [axis('H', 9, 0.1), axis('K', 9, 0.1), axis('L', 7, 0.1)];
  const shape = [7, 9, 9], N = 7 * 9 * 9;
  const volume = randomVolume(N, 3);
  const mask = Uint8Array.from({ length: N }, (_, i) => (i % 13 === 0 ? 1 : 0));
  const rand = random(5), variance = Float32Array.from({ length: N }, () => rand());
  for (const split of [1, 2]) {
    const plan = powderPlan(dims, hexCell, { bins: [0.05], split });
    for (const [name] of PRESETS) {
      const maps = indexMaps(group(name), dims);
      const got = powderAverage(volume, shape, plan, maps, mask, variance);
      const want = reference(volume, shape, plan, maps, mask, variance);
      for (let b = 0; b < plan.count; b++) {
        const label = `${name} split ${split} shell ${b}`;
        assert.ok(close(got.voxels[b], want.voxels[b]), `${label} voxels ${got.voxels[b]} vs ${want.voxels[b]}`);
        assert.ok(close(got.intensity[b], want.intensity[b]), `${label} I ${got.intensity[b]} vs ${want.intensity[b]}`);
        assert.ok(close(got.sigma[b], want.sigma[b]), `${label} σ ${got.sigma[b]} vs ${want.sigma[b]}`);
      }
    }
  }
});

test('orbits count with their multiplicity, not their measured members', () => {
  // Cubic a = 1 Å under mmm: the six voxels at |h| = 1 form three orbits of two.
  // (1,0,0) = 100 is measured but (-1,0,0) is not; the four others are 1.
  const dims = [axis('H', 3, 1), axis('K', 3, 1), axis('L', 3, 1)];
  const at = (h, k, l) => (h + 1) + 3 * ((k + 1) + 3 * (l + 1));
  const volume = new Float32Array(27).fill(NaN), variance = new Float32Array(27);
  volume[at(1, 0, 0)] = 100;
  variance[at(1, 0, 0)] = 4;
  for (const [h, k, l] of [[0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]) {
    volume[at(h, k, l)] = 1;
    variance[at(h, k, l)] = 1;
  }
  const plan = powderPlan(dims, cubic, { bins: [1], split: 1 });
  const shell = Math.floor(2 * Math.PI);
  const result = powderAverage(volume, [3, 3, 3], plan, indexMaps(group('mmm'), dims), null, variance);
  // The orbit of (1,0,0) stands for two voxels of the shell: (2·100 + 4·1) / 6,
  // not (100 + 4·1) / 5 as pooling the measured voxels would give.
  assert.ok(close(result.intensity[shell], 34));
  assert.equal(result.voxels[shell], 6);
  // Both members of an orbit carry its one mean: σ² = (2²·4 + 2²·½ + 2²·½) / 6².
  assert.ok(close(result.sigma[shell], Math.sqrt(20 / 36)));
  // Without symmetry, the five measured voxels are averaged.
  const plain = powderAverage(volume, [3, 3, 3], plan, [IDENTITY_MAP], null, variance);
  assert.ok(close(plain.intensity[shell], 20.8));
  assert.ok(close(plain.sigma[shell], Math.sqrt(8) / 5));
});

test('coverage gaps and masks change the precision, not the intensity', () => {
  const dims = [axis('H', 41, 0.1), axis('K', 41, 0.1), axis('L', 41, 0.1)];
  const G = reciprocalMetric(hexCell), q2 = (x) => 4 * Math.PI ** 2 * [0, 1, 2].reduce((s, i) => s + x[i] * [0, 1, 2].reduce((t, j) => t + G[i][j] * x[j], 0), 0);
  // A smooth isotropic intensity, and a constant.
  const { volume, shape } = fill(dims, (x) => 5 + 100 * Math.exp(-q2(x) / 8));
  const flat = fill(dims, () => 7).volume;
  const rand = random(11), holes = Uint8Array.from(volume, () => (rand() < 0.4 ? 1 : 0));
  const mask = Uint8Array.from(volume, () => (rand() < 0.1 ? 1 : 0));
  const gappy = Float32Array.from(volume, (v, i) => (holes[i] ? NaN : v));
  const flatGappy = Float32Array.from(flat, (v, i) => (holes[i] ? NaN : v));
  for (const split of [1, 2, 3]) {
    const plan = powderPlan(dims, hexCell, { split });
    const full = powderAverage(volume, shape, plan, [IDENTITY_MAP]);
    const partial = powderAverage(gappy, shape, plan, [IDENTITY_MAP], mask);
    const constant = powderAverage(flatGappy, shape, plan, [IDENTITY_MAP], mask);
    // Shells inside the grid's inscribed sphere (the nearest face, L = 2.05, is
    // 2π·2.05/c from the origin), away from the origin.
    const inner = 0.9 * 2 * Math.PI * 2.05 / 6.8;
    let checked = 0, kept = 0, all = 0, covered = 0;
    for (let b = 3; plan.edges[b + 1] < inner; b++, checked++) {
      const q = (plan.edges[b] + plan.edges[b + 1]) / 2, label = `split ${split} shell ${b}`;
      assert.ok(Math.abs(full.intensity[b] / (5 + 100 * Math.exp(-q * q / 8)) - 1) < 0.02, `${label}: ${full.intensity[b]}`);
      assert.ok(Math.abs(partial.intensity[b] / full.intensity[b] - 1) < 0.02, `${label}: ${partial.intensity[b]} vs ${full.intensity[b]}`);
      assert.ok(close(constant.intensity[b], 7, 1e-6), label);
      // A full shell is covered once. Voxel centres alias against thin shells
      // (±30% here); split voxels follow the shell volume.
      if (split > 1) assert.ok(Math.abs(full.coverage[b] - 1) < 0.05, `${label} coverage ${full.coverage[b]}`);
      covered += full.coverage[b];
      kept += partial.voxels[b];
      all += full.voxels[b];
    }
    assert.ok(checked >= 10);
    assert.ok(Math.abs(covered / checked - 1) < 0.03, `split ${split}: mean coverage ${covered / checked}`);
    // 60% of the voxels are left by the holes, and 90% of those by the mask.
    assert.ok(Math.abs(kept / all - 0.54) < 0.03, `split ${split}: ${kept / all}`);
  }
});

test('split voxels share each voxel between shells without losing any', () => {
  const dims = [axis('H', 9, 0.1), axis('K', 9, 0.1), axis('L', 7, 0.2)];
  const volume = randomVolume(7 * 9 * 9, 17), measured = volume.filter((v) => !Number.isNaN(v)).length;
  for (const split of [1, 2, 3, 4]) {
    const plan = powderPlan(dims, hexCell, { bins: [0.03], split });
    const { voxels } = powderAverage(volume, [7, 9, 9], plan, [IDENTITY_MAP]);
    assert.ok(Math.abs(voxels.reduce((a, b) => a + b) - measured) < 1e-9, `split ${split}`);
  }
  // A shell narrower than a voxel gets a share of it with split voxels only.
  const one = new Float32Array(7 * 9 * 9).fill(NaN);
  one[6 + 9 * (4 + 9 * 3)] = 1;
  const coarse = powderAverage(one, [7, 9, 9], powderPlan(dims, hexCell, { bins: [0.01], split: 1 }), [IDENTITY_MAP]);
  const fine = powderAverage(one, [7, 9, 9], powderPlan(dims, hexCell, { bins: [0.01], split: 4 }), [IDENTITY_MAP]);
  assert.equal(coarse.voxels.filter((x) => x > 0).length, 1);
  assert.ok(fine.voxels.filter((x) => x > 0).length > 5);
});

test('Q bins are read as Mantid Rebin parameters', () => {
  assert.equal(parseBins(''), null);
  assert.equal(parseBins('  '), null);
  assert.deepEqual(parseBins('0.05'), [0.05]);
  assert.deepEqual(parseBins('-0.01'), [-0.01]);
  assert.deepEqual(parseBins('0.5, 0.02, 3, 0.05, 10'), [0.5, 0.02, 3, 0.05, 10]);
  assert.deepEqual(parseBins('0.5 0.02 3'), [0.5, 0.02, 3]);
  assert.throws(() => parseBins('0.05 Å'), /numbers/);
  assert.throws(() => parseBins('0.5, 0.02'), /one step, or ranges/);
  assert.throws(() => parseBins('1, 0, 2'), /cannot be 0/);
  assert.throws(() => parseBins('0'), /cannot be 0/);
  assert.throws(() => parseBins('3, 0.1, 2'), /increase/);
  assert.throws(() => parseBins('-1, 0.1, 2'), /0 or above/);
});

test('shell edges follow Mantid Rebin', () => {
  const round = (edges) => Array.from(edges, (e) => +e.toFixed(9));
  // Uniform edges are counted from the start, so they do not drift.
  const tenth = shellEdges([0.1], 0, 30);
  assert.equal(tenth.length, 301);
  assert.ok(tenth.every((e, b) => Math.abs(e - b / 10) < 1e-12));
  // A range ends at its boundary with a last bin of 0.25 to 1.25 steps.
  assert.deepEqual(round(shellEdges([0.3], 0, 1)), [0, 0.3, 0.6, 0.9, 1]);
  assert.deepEqual(round(shellEdges([0.45], 0, 1)), [0, 0.45, 1]);
  // Logarithmic: each edge 1.1 times the previous, up to the boundary.
  const log = shellEdges([-0.1], 1, 2);
  assert.equal(log.at(-1), 2);
  for (let b = 1; b < log.length - 1; b++) assert.ok(Math.abs(log[b] / log[b - 1] - 1.1) < 1e-12, `edge ${b}`);
  assert.ok(log.at(-1) - log.at(-2) >= 0.25 * 0.1 * log.at(-2));
  // Ranges with their own steps, uniform then logarithmic.
  assert.deepEqual(round(shellEdges([0.5, 0.1, 1, 0.25, 2])), [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2]);
  const mixed = shellEdges([0.5, 0.25, 1, -0.5, 4]);
  assert.deepEqual(round(mixed), [0.5, 0.75, 1, 1.5, 2.25, 3.375, 4]);
  assert.throws(() => shellEdges([-0.01], 0, 5), /start above 0/);
  assert.throws(() => shellEdges([0.1], 2, 1), /above Q min/);
});

test('uneven shells match the direct implementation and keep every voxel', () => {
  const dims = [axis('H', 9, 0.1), axis('K', 9, 0.1), axis('L', 7, 0.1)];
  const shape = [7, 9, 9], N = 7 * 9 * 9;
  const volume = randomVolume(N, 23);
  const rand = random(29), variance = Float32Array.from({ length: N }, () => rand());
  for (const bins of [[-0.04], [0.2, 0.03, 0.8, 0.2, 2.5], [0, 0.5, 1, -0.1, 3]]) {
    for (const name of ['1', '6/mmm', 'm-3m']) {
      const plan = powderPlan(dims, hexCell, { bins, split: 2 });
      const maps = indexMaps(group(name), dims);
      const got = powderAverage(volume, shape, plan, maps, null, variance);
      const want = reference(volume, shape, plan, maps, null, variance);
      for (let b = 0; b < plan.count; b++) {
        const label = `${bins} ${name} shell ${b}`;
        assert.ok(close(got.voxels[b], want.voxels[b]), `${label} voxels ${got.voxels[b]} vs ${want.voxels[b]}`);
        assert.ok(close(got.intensity[b], want.intensity[b]), `${label} I`);
        assert.ok(close(got.sigma[b], want.sigma[b]), `${label} σ`);
      }
    }
  }
  // Without symmetry, the shells from Q min to Q max hold exactly the voxels whose centres lie there.
  const plan = powderPlan(dims, hexCell, { bins: [0.2, 0.07, 0.9, -0.2, 2], split: 1 });
  const { voxels } = powderAverage(volume, shape, plan, [IDENTITY_MAP]);
  const qOf = (x) => Math.sqrt(x.reduce((s, xi, i) => s + xi * x.reduce((t, xj, j) => t + plan.G[3 * i + j] * xj, 0), 0));
  let inside = 0;
  for (let l = 0; l < 7; l++) for (let k = 0; k < 9; k++) for (let h = 0; h < 9; h++) {
    const q = qOf([plan.c0[0] + h * 0.1, plan.c0[1] + k * 0.1, plan.c0[2] + l * 0.1]);
    if (!Number.isNaN(volume[h + 9 * (k + 9 * l)]) && q >= 0.2 && q < 2) inside++;
  }
  assert.ok(inside > 50);
  assert.ok(Math.abs(voxels.reduce((a, b) => a + b) - inside) < 1e-9);
  // Coverage uses each shell's own volume: shells of 0.1 and 0.2 Å⁻¹ inside a full grid are covered once.
  const big = [axis('H', 41, 0.1), axis('K', 41, 0.1), axis('L', 41, 0.1)];
  const flat = powderAverage(new Float32Array(41 ** 3).fill(1), [41, 41, 41], powderPlan(big, hexCell, { bins: [0.3, 0.1, 0.8, 0.2, 1.6] }), [IDENTITY_MAP]);
  assert.ok(flat.coverage.every((f) => Math.abs(f - 1) < 0.05), Array.from(flat.coverage).join(' '));
});

test('shells, metric and refusals', () => {
  const dims = [axis('H', 101, 0.1), axis('K', 101, 0.1), axis('L', 101, 0.1)];
  const plan = powderPlan(dims, hexCell);
  // Shortest step: 0.1 r.l.u. along L, 2π/6.8 × 0.1 = 0.0924 Å⁻¹.
  assert.deepEqual(plan.bins, [0.09]);
  assert.ok(Math.abs(plan.step - 2 * Math.PI * 0.1 / 6.8) < 1e-12);
  // Farthest corner (5.05, 5.05, ±5.05): a*²(h² + k² + hk) + c*² l², with 2π.
  const astar2 = 4 / (3 * 4.2 ** 2), top = 2 * Math.PI * Math.sqrt(astar2 * 3 * 5.05 ** 2 + 5.05 ** 2 / 6.8 ** 2);
  assert.ok(Math.abs(plan.top - top) < 1e-9);
  // Shells of 0.09 from 0, the last one ending at the farthest corner.
  assert.equal(plan.edges[0], 0);
  assert.equal(plan.edges.at(-1), plan.top);
  assert.ok(plan.edges.slice(0, -1).every((e, b) => Math.abs(e - 0.09 * b) < 1e-12));
  // One voxel is (2π)³ Δh Δk Δl / V_cell.
  const cellVolume = 4.2 * 4.2 * 6.8 * Math.sin(Math.PI / 3);
  assert.ok(Math.abs(plan.voxel * cellVolume / (8 * Math.PI ** 3 * 1e-3) - 1) < 1e-9);
  // Q min and Q max bound a single step; ranges in the bins override them.
  assert.deepEqual(Array.from(powderPlan(dims, hexCell, { bins: [0.2], qmin: 1, qmax: 1.62 }).edges, (e) => +e.toFixed(9)), [1, 1.2, 1.4, 1.62]);
  assert.deepEqual(Array.from(powderPlan(dims, hexCell, { bins: [1, 0.5, 2], qmin: 0, qmax: 9 }).edges), [1, 1.5, 2]);
  // Logarithmic bins start at the shortest step by default.
  assert.equal(powderPlan(dims, hexCell, { bins: [-0.05] }).edges[0], 0.092);
  assert.equal(powderPlan(dims, hexCell, { edges: [0, 1, 3] }).count, 2);

  // Oblique axes: [H,H,0] at H = 1 is (110), |Q| = 2π√3 a* (γ* = 60°), d = a/2.
  const hh = { ...dims[0], label: '[H,H,0]', basis: parseBasis('[H,H,0]') };
  const { G } = qMetric([hh, dims[1], dims[2]], hexCell);
  assert.ok(Math.abs(Math.sqrt(G[0]) - 2 * Math.PI / (4.2 / 2)) < 1e-12);

  // Q frames in Å⁻¹ are Cartesian.
  const q = (label) => ({ label, units: 'Angstrom^-1', frame: 'QSample', basis: null, edges: Float64Array.of(-1, 0, 1) });
  assert.deepEqual(qMetric([q('Qx'), q('Qy'), q('Qz')], null), { G: [1, 0, 0, 0, 1, 0, 0, 0, 1], frame: 'Q' });
  for (const units of ['Å⁻¹', '1/Angstrom', 'A^-1', 'inverse angstroms']) assert.equal(qMetric([q('a'), q('b'), { ...q('c'), units }], null).frame, 'Q', units);

  assert.throws(() => powderPlan(dims, null), /unit cell/);
  const energy = { label: 'DeltaE', units: 'meV', frame: 'General Frame', basis: null, edges: dims[2].edges };
  assert.throws(() => powderPlan([dims[0], dims[1], energy], hexCell), /H, K, L axes/);
  const uneven = { ...dims[2], edges: Float64Array.from(dims[2].edges, (e, i) => e + (i === 3 ? 0.03 : 0)) };
  assert.throws(() => powderPlan([dims[0], dims[1], uneven], hexCell), /non-uniform/);
  assert.throws(() => powderPlan(dims, hexCell, { bins: [1e-5] }), /shells/);
  assert.throws(() => powderPlan(dims, hexCell, { split: 5 }), /sub-cells/);
});
