import assert from 'node:assert/strict';
import test from 'node:test';

import { EDGE, edgeMask, maskStats, OUTLIER, outlierMask } from '../js/mask.js';
import { averageSlab } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps } from '../js/symmetry.js';

function randomVolume(shape, seed, holes = 0.1) {
  let s = seed;
  const rand = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648);
  return Float32Array.from({ length: shape[0] * shape[1] * shape[2] }, () => (rand() < holes ? NaN : 10 + rand()));
}

test('edgeMask flags measured voxels within the box radius of unmeasured ones', () => {
  const shape = [6, 7, 9], volume = randomVolume(shape, 3, 0.04);
  for (const radius of [1, 2]) {
    const mask = edgeMask(volume, shape, radius);
    for (let a = 0; a < shape[0]; a++) for (let b = 0; b < shape[1]; b++) for (let c = 0; c < shape[2]; c++) {
      const i = (a * shape[1] + b) * shape[2] + c;
      let near = false;
      for (let da = -radius; da <= radius; da++) for (let db = -radius; db <= radius; db++) for (let dc = -radius; dc <= radius; dc++) {
        const [x, y, z] = [a + da, b + db, c + dc];
        if (x < 0 || y < 0 || z < 0 || x >= shape[0] || y >= shape[1] || z >= shape[2]) continue;
        if (Number.isNaN(volume[(x * shape[1] + y) * shape[2] + z])) near = true;
      }
      assert.equal(mask[i], near && !Number.isNaN(volume[i]) ? EDGE : 0, `voxel ${a},${b},${c} radius ${radius}`);
    }
  }
  assert.ok(edgeMask(volume, shape, 0).every((m) => m === 0));
});

test('outlierMask rejects a spike that disagrees with its equivalents', () => {
  // Storage (L, K, H) = (3, 5, 5); H and K share a width, so 6/m maps the grid.
  const shape = [3, 5, 5];
  const axis = (n, w, name) => ({
    label: name, basis: { vec: [...'HKL'].map((c) => (c === name ? 1 : 0)) },
    edges: Array.from({ length: n + 1 }, (_, i) => (i - n / 2) * w),
  });
  const maps = indexMaps(closeGroup(parseOps('h+k,-h,l; -h,-k,-l')), [axis(5, 0.1, 'H'), axis(5, 0.1, 'K'), axis(3, 0.2, 'L')]);
  const volume = randomVolume(shape, 7, 0);
  const clean = outlierMask(volume, shape, maps, 5);
  assert.equal(clean.filter((m) => m).length, 0);

  const spike = (2 * 5 + 1) * 5 + 3; // L index 2, K index 1, H index 3
  volume[spike] = 1000;
  const mask = outlierMask(volume, shape, maps, 5);
  assert.deepEqual([...mask.keys()].filter((i) => mask[i]), [spike]);
  assert.equal(mask[spike], OUTLIER);
  assert.deepEqual(maskStats(volume, mask), { measured: 75, edge: 0, outlier: 1 });

  // Masked voxels drop out of the averages, and invert shows only them.
  const kept = averageSlab(volume, shape, 2, [2], maps, mask);
  const removed = averageSlab(volume, shape, 2, [2], maps, mask, true);
  assert.ok(Math.max(...kept.values.filter(Number.isFinite)) < 20);
  assert.equal(removed.values.filter(Number.isFinite).length > 0, true);
  assert.ok(removed.values.filter(Number.isFinite).every((v) => v === 1000));
});
