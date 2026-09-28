import assert from 'node:assert/strict';
import test from 'node:test';

import { binVolume, coarseGrid, orbitMean, surfaceNets } from '../js/iso.js';
import { closeGroup, indexMaps, parseOps } from '../js/symmetry.js';

const axis = (n, w, name) => ({
  label: name, basis: { vec: [...'HKL'].map((c) => (c === name ? 1 : 0)), letter: name },
  edges: Array.from({ length: n + 1 }, (_, i) => (i - n / 2) * w),
});

test('coarse blocks stay centered on the origin', () => {
  const grid = coarseGrid([axis(401, 0.1, 'H'), axis(400, 0.1, 'K'), axis(11, 0.25, 'L')], 100);
  assert.equal(grid.factor, 5);
  assert.deepEqual(grid.shape, [81, 80, 3]);
  // Odd axis: a block centered at 0; even axis: a block edge at 0.
  const [h, k] = grid.dims.map((d) => d.edges);
  assert.ok(Math.abs(h[40] + h[41]) < 1e-9);
  assert.ok(Math.abs(k[40]) < 1e-9);
  assert.ok(Math.abs(h[0] + h[81]) < 1e-9 && Math.abs(k[0] + k[80]) < 1e-9);
  // The middle fine bin lands in the middle block.
  assert.equal(grid.block[0][200], 40);
});

test('binVolume sums valid voxels per block', () => {
  const dims = [axis(5, 0.1, 'H'), axis(3, 0.1, 'K'), axis(3, 0.1, 'L')];
  const grid = coarseGrid(dims, 2);
  assert.equal(grid.factor, 3);
  const volume = Float32Array.from({ length: 45 }, (_, i) => (i % 4 ? 1 : NaN));
  const { sums, counts } = binVolume(volume, [3, 3, 5], grid);
  assert.equal(counts.reduce((a, b) => a + b, 0), 45 - 12);
  assert.equal(sums.reduce((a, b) => a + b, 0), 45 - 12);
});

test('orbitMean pools the sums and counts of each orbit', () => {
  const shape = [5, 5, 3];
  const dims = [axis(5, 0.1, 'H'), axis(5, 0.1, 'K'), axis(3, 0.2, 'L')];
  const ops = closeGroup(parseOps('h+k,-h,l; -h,-k,-l'));
  const maps = indexMaps(ops, dims);
  const sums = Float64Array.from({ length: 75 }, (_, i) => (i * 7919) % 13);
  const counts = Uint32Array.from({ length: 75 }, (_, i) => (i * 31) % 3);
  const mean = orbitMean({ sums, counts }, shape, maps);
  for (let p = 0; p < 75; p++) {
    const i = [p % 5, Math.floor(p / 5) % 5, Math.floor(p / 25)];
    const orbit = new Set();
    for (const { M, t } of maps) {
      const q = [0, 1, 2].map((r) => M[3 * r] * i[0] + M[3 * r + 1] * i[1] + M[3 * r + 2] * i[2] + t[r]);
      if (q.every((x, d) => x >= 0 && x < shape[d])) orbit.add(q[0] + 5 * q[1] + 25 * q[2]);
    }
    let s = 0, c = 0;
    for (const q of orbit) { s += sums[q]; c += counts[q]; }
    if (c) assert.ok(Math.abs(mean[p] - s / c) < 1e-6, `voxel ${p}`);
    else assert.ok(Number.isNaN(mean[p]));
  }
});

test('surfaceNets draws a closed sphere', () => {
  const n = 24, r = 8, c = (n - 1) / 2;
  const field = new Float32Array(n ** 3);
  for (let k = 0, i = 0; k < n; k++) for (let j = 0; j < n; j++) for (let h = 0; h < n; h++, i++) {
    field[i] = r - Math.hypot(h - c, j - c, k - c);
  }
  const { positions, indices } = surfaceNets(field, [n, n, n], 0);
  assert.ok(indices.length > 1000);
  for (let v = 0; v < positions.length; v += 3) {
    const d = Math.hypot(positions[v] - c, positions[v + 1] - c, positions[v + 2] - c);
    assert.ok(Math.abs(d - r) < 0.5, `vertex at radius ${d}`);
  }
  // Closed and consistently oriented: every directed edge appears once, reversed once.
  const edges = new Map();
  for (let t = 0; t < indices.length; t += 3) {
    for (let e = 0; e < 3; e++) {
      const a = indices[t + e], b = indices[t + (e + 1) % 3];
      edges.set(`${a},${b}`, (edges.get(`${a},${b}`) ?? 0) + 1);
    }
  }
  for (const [key, count] of edges) {
    const [a, b] = key.split(',');
    assert.equal(count, 1, key);
    assert.equal(edges.get(`${b},${a}`), 1, `reverse of ${key}`);
  }
  // No data counts as outside.
  field.fill(NaN, 0, n * n * 12);
  assert.ok(surfaceNets(field, [n, n, n], 0).indices.length > 0);
});
