import assert from 'node:assert/strict';
import { mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { exportPlan, nebulaUB, sigmaNote, symmetrizeForExport, writeNebulaFile } from '../js/export.js';
import { reciprocalMetric } from '../js/nexus.js';
import { averageSlab, IDENTITY_MAP } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps, PRESETS } from '../js/symmetry.js';

const group = (name) => closeGroup(parseOps(PRESETS.find(([n]) => n === name)[1]));
// An HKL axis with n bins of width w starting at `start`.
const axis = (letter, n, w, start) => ({
  label: letter, basis: { vec: [...'HKL'].map((c) => (c === letter ? 1 : 0)), letter },
  edges: Float64Array.from({ length: n + 1 }, (_, i) => start + i * w),
});
const hexCell = { a: 4.2, b: 4.2, c: 6.8, alpha: 90, beta: 90, gamma: 120 };

function randomVolume(n, seed = 1) {
  let s = seed;
  const rand = () => ((s = (s * 16807) % 2147483647) / 2147483647);
  return Float32Array.from({ length: n }, () => (rand() < 0.3 ? NaN : rand() * 10));
}

test('plan: symmetric grids are kept, others padded about 0', () => {
  const dims = [axis('H', 7, 0.1, -0.35), axis('K', 7, 0.1, -0.35), axis('L', 5, 0.2, -0.5)];
  const plan = exportPlan(dims, hexCell);
  assert.deepEqual(plan.shape, [7, 7, 5]);
  assert.equal(plan.padded, false);
  // Bin centred at 0, asymmetric range: odd grid with 0 in the middle.
  const odd = exportPlan([axis('H', 6, 0.1, -0.15), dims[1], dims[2]], hexCell);
  assert.deepEqual([odd.size[0], odd.lo[0]], [9, -3]);
  assert.ok(odd.centers[0].every((c, i) => Math.abs(c - (i - 4) * 0.1) < 1e-12));
  // Display order K, L, H is reordered to H, K, L.
  const permuted = exportPlan([dims[1], dims[2], dims[0]], hexCell);
  assert.deepEqual(permuted.order, [2, 0, 1]);
});

test('plan: NEBULA3D\'s origin, the centre of bin n//2, is Q = 0 on every axis', () => {
  // NEBULA3D's FFT takes index n//2 of each axis as Q = 0 (ifftshift).
  const K = axis('K', 7, 0.1, -0.35), L = axis('L', 5, 0.2, -0.5);
  for (const [n, start] of [[7, -0.35], [6, -0.15], [9, -0.25], [4, -0.05], [11, -0.85]]) {
    const plan = exportPlan([axis('H', n, 0.1, start), K, L], hexCell);
    plan.centers.forEach((c, a) => {
      assert.equal(c.length % 2, 1, `axis ${a} of H from ${start}: odd`);
      assert.ok(Math.abs(c[c.length >> 1]) < 1e-12, `axis ${a} of H from ${start}: centre ${c[c.length >> 1]}`);
    });
  }
  // A bin edge at 0 (an even grid split there) would put +w/2 at n//2: refused, with the binning that works.
  // The reviewer's case: H from -2 to 2 in 40 bins of 0.1.
  for (const [n, start] of [[6, -0.2], [40, -2], [5, 0]]) {
    assert.throws(() => exportPlan([axis('H', n, 0.1, start), K, L], hexCell), /H has no bin centred at 0.*bin edge lies at 0/);
  }
  assert.throws(() => exportPlan([axis('H', 40, 0.1, -2), K, L], hexCell), /-2\.05 to 2\.05 in 41 bins of 0\.1 \(Mantid binning "-2\.05,0\.1,2\.05"\)/);
});

test('plan: data that cannot be exported is refused', () => {
  const dims = [axis('H', 7, 0.1, -0.35), axis('K', 7, 0.1, -0.35), axis('L', 5, 0.2, -0.5)];
  assert.throws(() => exportPlan(dims, null), /unit cell/);
  const projected = { ...dims[0], label: '[H,H,0]', basis: { vec: [1, 1, 0], letter: 'H' } };
  assert.throws(() => exportPlan([projected, dims[1], dims[2]], hexCell), /plain H, K and L/);
  const uneven = { ...dims[0], edges: Float64Array.of(-0.35, -0.25, -0.1, 0.05, 0.15, 0.25, 0.35, 0.45) };
  assert.throws(() => exportPlan([uneven, dims[1], dims[2]], hexCell), /non-uniform/);
  assert.throws(() => exportPlan([axis('H', 6, 0.1, -0.23), dims[1], dims[2]], hexCell), /no bin centred at 0.*±0\.05, ±0\.15/);
  assert.throws(() => exportPlan([axis('H', 4, 0.1, 1.05), dims[1], dims[2]], hexCell), /too far/);
});

test('UB for NEBULA3D carries 2π and the reciprocal metric', () => {
  const ub = [0.1, 0.2, 0.05, -0.2, 0.1, 0, 0.03, 0, 0.14];
  assert.deepEqual(nebulaUB({ ...hexCell, ub }), ub.map((x) => 2 * Math.PI * x));
  // Without a UB: (UB)ᵀUB = (2π)² G*.
  const B = nebulaUB(hexCell), G = reciprocalMetric(hexCell);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const btb = B[i] * B[j] + B[3 + i] * B[3 + j] + B[6 + i] * B[6 + j];
      assert.ok(Math.abs(btb - 4 * Math.PI ** 2 * G[i][j]) < 1e-9, `${i}${j}`);
    }
  }
});

test('exported planes equal one-bin symmetrized slices', () => {
  const dims = [axis('H', 9, 0.1, -0.45), axis('K', 9, 0.1, -0.45), axis('L', 7, 0.2, -0.7)];
  const shape = [7, 9, 9], volume = randomVolume(7 * 9 * 9);
  const mask = Uint8Array.from(volume, (_, i) => (i % 11 === 0 ? 1 : 0));
  for (const name of ['1', '-1', 'mmm', '6/mmm']) {
    const plan = exportPlan(dims, hexCell);
    const maps = indexMaps(group(name), plan.paddedDims);
    const { data, valid } = symmetrizeForExport(volume, shape, plan, maps, { mask });
    for (let l = 0; l < 7; l++) {
      const slab = averageSlab(volume, shape, 2, [l], maps, mask);
      for (let k = 0; k < 9; k++) {
        for (let h = 0; h < 9; h++) {
          const v = slab.values[k * 9 + h], o = (h * 9 + k) * 7 + l;
          if (Number.isNaN(v)) assert.equal(valid[o], 0, `${name} ${h} ${k} ${l}`);
          else assert.ok(valid[o] === 1 && Math.abs(data[o] - v) < 1e-4 * Math.max(1, Math.abs(v)), `${name} ${h} ${k} ${l}`);
        }
      }
    }
  }
});

test('σ of each exported voxel is that of its orbit\'s mean, √(Σσ²)/n over its distinct valid voxels', () => {
  // H padded from -0.25 … 0.45 to ±0.45, so orbits meet the padding and the mirror planes.
  const dims = [axis('H', 8, 0.1, -0.35), axis('K', 9, 0.1, -0.45), axis('L', 7, 0.2, -0.7)];
  const shape = [7, 9, 8], volume = randomVolume(7 * 9 * 8, 3);
  const variance = Float32Array.from(randomVolume(volume.length, 5), (v, i) => (i % 13 === 0 ? NaN : Number.isNaN(v) ? 0.5 : v / 10));
  const mask = Uint8Array.from(volume, (_, i) => (i % 7 === 0 ? 1 : 0));
  const plan = exportPlan(dims, hexCell), [P0, P1, P2] = plan.size;
  for (const name of ['1', '-1', 'mmm', '6/mmm']) {
    const maps = indexMaps(group(name), plan.paddedDims);
    const { data, valid, sigma, stats } = symmetrizeForExport(volume, shape, plan, maps, { mask, variance });
    let defaulted = 0;
    for (let p2 = 0; p2 < P2; p2++) for (let p1 = 0; p1 < P1; p1++) for (let p0 = 0; p0 < P0; p0++) {
      // The orbit's distinct valid source voxels, by set.
      const members = new Set();
      for (const { M, t } of maps) {
        const q = [0, 1, 2].map((i) => M[3 * i] * p0 + M[3 * i + 1] * p1 + M[3 * i + 2] * p2 + t[i]);
        const s = q.map((qi, i) => qi + plan.lo[i]);
        if (q.some((qi, i) => qi < 0 || qi >= plan.size[i]) || s.some((si, i) => si < 0 || si >= [8, 9, 7][i])) continue;
        const f = s[0] + 8 * (s[1] + 9 * s[2]);
        if (!Number.isNaN(volume[f]) && !mask[f]) members.add(f);
      }
      const o = p0 * 9 * 7 + p1 * 7 + p2;
      if (!members.size) {
        assert.equal(valid[o], 0);
        assert.equal(sigma[o], 0);
        continue;
      }
      const known = [...members].filter((f) => variance[f] >= 0);
      const mean = [...members].reduce((a, f) => a + volume[f], 0) / members.size;
      const expected = known.length ? Math.sqrt(known.reduce((a, f) => a + variance[f], 0)) / known.length : Math.sqrt(Math.abs(mean));
      if (!known.length) defaulted++;
      assert.ok(valid[o] === 1 && Math.abs(data[o] - mean) < 1e-4 * Math.max(1, Math.abs(mean)), `${name} value ${p0} ${p1} ${p2}`);
      assert.ok(Math.abs(sigma[o] - expected) < 1e-5 * Math.max(1, expected), `${name} σ ${p0} ${p1} ${p2}: ${sigma[o]} vs ${expected}`);
    }
    assert.equal(stats.sigmaDefaulted, defaulted, name);
  }
  // Without uncertainties there is no sigma, and the file says NEBULA3D takes √|data|.
  const plain = symmetrizeForExport(volume, shape, plan, [IDENTITY_MAP], { mask });
  assert.equal(plain.sigma, null);
  assert.match(sigmaNote(plain, undefined), /none written.*sqrt\(\|data\|\)/);
});

test('symmetry fills the padding from measured equivalents', () => {
  // H covers -0.2 … 0.4 (bins centred at -0.2 … 0.4): padded to -0.4 … 0.4.
  const dims = [axis('H', 7, 0.1, -0.25), axis('K', 3, 0.1, -0.15), axis('L', 3, 0.1, -0.15)];
  const shape = [3, 3, 7], volume = Float32Array.from({ length: 63 }, (_, i) => i);
  const plan = exportPlan(dims, hexCell);
  assert.deepEqual(plan.shape, [9, 3, 3]);
  const at = (h, k, l) => (h * 3 + k) * 3 + l;
  // Without symmetry the padding is invalid; with inversion it mirrors +H.
  const plain = symmetrizeForExport(volume, shape, plan, [IDENTITY_MAP]);
  assert.equal(plain.valid[at(0, 1, 1)], 0);
  assert.equal(plain.data[at(8, 1, 1)], volume[6 + 7 * (1 + 3 * 1)]);
  const inv = symmetrizeForExport(volume, shape, plan, indexMaps(group('-1'), plan.paddedDims));
  // (-0.4, 0, 0) pairs only with (0.4, 0, 0), so it takes that voxel's value.
  assert.equal(inv.valid[at(0, 1, 1)], 1);
  assert.equal(inv.data[at(0, 1, 1)], volume[6 + 7 * (1 + 3 * 1)]);
  assert.equal(inv.stats.valid, 81);
});

test('the file has NEBULA3D\'s layout', async () => {
  await h5wasm.ready;
  const dir = mkdtempSync(join(tmpdir(), 'nxv-export-'));
  try {
    const dims = [axis('H', 5, 0.1, -0.25), axis('K', 5, 0.1, -0.25), axis('L', 3, 0.2, -0.3)];
    const shape = [3, 5, 5], volume = randomVolume(75, 7);
    const plan = exportPlan(dims, { ...hexCell, ub: [0.27, 0.14, 0, 0, 0.24, 0, 0, 0, 0.147] });
    const variance = Float32Array.from(volume, (v) => Math.abs(v) / 4);
    const result = symmetrizeForExport(volume, shape, plan, [IDENTITY_MAP], { variance });
    const path = join(dir, 'out.nxs');
    writeNebulaFile(h5wasm, path, plan, result, { source_file: 'x.nxs', symmetry: '1' });
    const f = new h5wasm.File(path, 'r');
    const entry = f.get('entry');
    assert.deepEqual(entry.get('data').shape, [5, 5, 3]);
    assert.deepEqual(entry.get('mask').shape, [5, 5, 3]);
    assert.deepEqual(Array.from(entry.get('h_axis').value, (x) => +x.toFixed(6)), [-0.2, -0.1, 0, 0.1, 0.2]);
    assert.deepEqual(entry.get('l_axis').value.length, 3);
    assert.ok(Math.abs(entry.get('ub_matrix').value[0] - 2 * Math.PI * 0.27) < 1e-12);
    assert.equal(entry.attrs.symmetry.value, '1');
    assert.equal(entry.attrs.instrument.value, '');
    // (H, K, L) order: data[h][k][l] is the source voxel at storage [l][k][h].
    const data = entry.get('data').value, mask = entry.get('mask').value;
    // σ as NEBULA3D reads it: /entry/sigma, the data's shape, 0 where invalid.
    assert.deepEqual(entry.get('sigma').shape, [5, 5, 3]);
    const sigma = entry.get('sigma').value;
    const v = volume[3 + 5 * (1 + 5 * 2)];
    const o = (3 * 5 + 1) * 3 + 2;
    if (Number.isNaN(v)) assert.ok(mask[o] === 0 && sigma[o] === 0);
    else assert.ok(mask[o] === 1 && Math.abs(data[o] - v) < 1e-6 && Math.abs(sigma[o] - Math.sqrt(Math.abs(v) / 4)) < 1e-6);
    assert.ok(mask.every((m, i) => m === 1 || sigma[i] === 0));
    f.close();
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
