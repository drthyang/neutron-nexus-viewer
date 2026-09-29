import assert from 'node:assert/strict';
import test from 'node:test';

import { parseBasis, reciprocalMetric } from '../js/nexus.js';
import { powderAverage, powderPlan, qMetric } from '../js/powder.js';
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
  const n = [shape[2], shape[1], shape[0]], { G, c0, w, dq, count, split } = plan;
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
        const shell = Math.floor(qOf(x) / dq);
        if (shell < count) shells.set(shell, (shells.get(shell) ?? 0) + 1 / split ** 3);
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
    const plan = powderPlan(dims, hexCell, { dq: 0.05, split });
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
  const plan = powderPlan(dims, cubic, { dq: 1, split: 1 });
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
    for (let b = 3; (b + 1) * plan.dq < inner; b++, checked++) {
      const q = (b + 0.5) * plan.dq, label = `split ${split} shell ${b}`;
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
    const plan = powderPlan(dims, hexCell, { dq: 0.03, split });
    const { voxels } = powderAverage(volume, [7, 9, 9], plan, [IDENTITY_MAP]);
    assert.ok(Math.abs(voxels.reduce((a, b) => a + b) - measured) < 1e-9, `split ${split}`);
  }
  // A shell narrower than a voxel gets a share of it with split voxels only.
  const one = new Float32Array(7 * 9 * 9).fill(NaN);
  one[6 + 9 * (4 + 9 * 3)] = 1;
  const coarse = powderAverage(one, [7, 9, 9], powderPlan(dims, hexCell, { dq: 0.01, split: 1 }), [IDENTITY_MAP]);
  const fine = powderAverage(one, [7, 9, 9], powderPlan(dims, hexCell, { dq: 0.01, split: 4 }), [IDENTITY_MAP]);
  assert.equal(coarse.voxels.filter((x) => x > 0).length, 1);
  assert.ok(fine.voxels.filter((x) => x > 0).length > 5);
});

test('shells, metric and refusals', () => {
  const dims = [axis('H', 101, 0.1), axis('K', 101, 0.1), axis('L', 101, 0.1)];
  const plan = powderPlan(dims, hexCell);
  // Shortest step: 0.1 r.l.u. along L, 2π/6.8 × 0.1 = 0.0924 Å⁻¹.
  assert.equal(plan.dq, 0.09);
  assert.ok(Math.abs(plan.step - 2 * Math.PI * 0.1 / 6.8) < 1e-12);
  // Farthest corner (5.05, 5.05, ±5.05): a*²(h² + k² + hk) + c*² l², with 2π.
  const astar2 = 4 / (3 * 4.2 ** 2), top = 2 * Math.PI * Math.sqrt(astar2 * 3 * 5.05 ** 2 + 5.05 ** 2 / 6.8 ** 2);
  assert.ok(Math.abs(plan.top - top) < 1e-9);
  assert.equal(plan.count, Math.ceil(top / 0.09));
  // One voxel is (2π)³ Δh Δk Δl / V_cell.
  const cellVolume = 4.2 * 4.2 * 6.8 * Math.sin(Math.PI / 3);
  assert.ok(Math.abs(plan.voxel * cellVolume / (8 * Math.PI ** 3 * 1e-3) - 1) < 1e-9);
  assert.deepEqual(powderPlan(dims, hexCell, { dq: 0.2, qmax: 3.05 }).count, 16);

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
  assert.throws(() => powderPlan(dims, hexCell, { dq: 1e-5 }), /shells/);
  assert.throws(() => powderPlan(dims, hexCell, { split: 5 }), /sub-cells/);
});
