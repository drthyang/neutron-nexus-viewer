import assert from 'node:assert/strict';
import { mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { exportPlan, nebulaUB, symmetrizeForExport, writeNebulaFile } from '../js/export.js';
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
  // Bin edge at 0: even grid split at 0.
  const even = exportPlan([axis('H', 6, 0.1, -0.2), dims[1], dims[2]], hexCell);
  assert.deepEqual([even.size[0], even.lo[0]], [8, -2]);
  assert.ok(Math.abs(even.centers[0][4] - 0.05) < 1e-12 && Math.abs(even.centers[0][3] + 0.05) < 1e-12);
  // Display order K, L, H is reordered to H, K, L.
  const permuted = exportPlan([dims[1], dims[2], dims[0]], hexCell);
  assert.deepEqual(permuted.order, [2, 0, 1]);
});

test('plan: data that cannot be exported is refused', () => {
  const dims = [axis('H', 7, 0.1, -0.35), axis('K', 7, 0.1, -0.35), axis('L', 5, 0.2, -0.5)];
  assert.throws(() => exportPlan(dims, null), /unit cell/);
  const projected = { ...dims[0], label: '[H,H,0]', basis: { vec: [1, 1, 0], letter: 'H' } };
  assert.throws(() => exportPlan([projected, dims[1], dims[2]], hexCell), /plain H, K and L/);
  const uneven = { ...dims[0], edges: Float64Array.of(-0.35, -0.25, -0.1, 0.05, 0.15, 0.25, 0.35, 0.45) };
  assert.throws(() => exportPlan([uneven, dims[1], dims[2]], hexCell), /non-uniform/);
  assert.throws(() => exportPlan([axis('H', 6, 0.1, -0.23), dims[1], dims[2]], hexCell), /neither/);
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
    const { data, valid } = symmetrizeForExport(volume, shape, plan, maps, mask);
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
    const result = symmetrizeForExport(volume, shape, plan, [IDENTITY_MAP]);
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
    const v = volume[3 + 5 * (1 + 5 * 2)];
    const o = (3 * 5 + 1) * 3 + 2;
    if (Number.isNaN(v)) assert.equal(mask[o], 0);
    else assert.ok(mask[o] === 1 && Math.abs(data[o] - v) < 1e-6);
    f.close();
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
