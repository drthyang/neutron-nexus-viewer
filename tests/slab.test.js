import assert from 'node:assert/strict';
import test from 'node:test';

import { averageSlab, IDENTITY_MAP, selectBins } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps, PRESETS } from '../js/symmetry.js';

// Direct statement of the semantics: pool the set of distinct in-range voxels
// g(p) over all operations g and slab voxels p under each output pixel.
function reference(volume, shape, fixed, ids, maps) {
  const n = [shape[2], shape[1], shape[0]];
  const [x, y] = [0, 1, 2].filter((d) => d !== fixed);
  const values = [], counts = [];
  for (let row = 0; row < n[y]; row++) {
    for (let col = 0; col < n[x]; col++) {
      const pooled = new Set();
      for (const j of ids) {
        const p = [0, 0, 0];
        p[x] = col; p[y] = row; p[fixed] = j;
        for (const { M, t } of maps) {
          const q = [0, 1, 2].map((i) => M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + t[i]);
          if (q.every((qi, i) => qi >= 0 && qi < n[i])) pooled.add(q[0] + n[0] * (q[1] + n[1] * q[2]));
        }
      }
      let total = 0, count = 0;
      for (const f of pooled) if (Number.isFinite(volume[f])) { total += volume[f]; count++; }
      counts.push(count);
      values.push(count ? total / count : NaN);
    }
  }
  return { values, counts };
}

function randomVolume(shape, seed) {
  let s = seed;
  const rand = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648);
  return Float32Array.from({ length: shape[0] * shape[1] * shape[2] }, () => (rand() < 0.15 ? NaN : rand() * 100 - 20));
}

// Display dims with origin-centered bins; H and K share a width.
function dims(nh, nk, nl) {
  const axis = (n, w, name) => ({
    label: name, basis: { vec: [...'HKL'].map((c) => (c === name ? 1 : 0)), letter: name },
    edges: Array.from({ length: n + 1 }, (_, i) => (i - n / 2) * w),
  });
  return [axis(nh, 0.1, 'H'), axis(nk, 0.1, 'K'), axis(nl, 0.25, 'L')];
}

const group = (name) => closeGroup(parseOps(PRESETS.find(([n]) => n === name)[1]));

function check(volume, shape, fixed, ids, maps, label) {
  const got = averageSlab(volume, shape, fixed, ids, maps);
  const want = reference(volume, shape, fixed, ids, maps);
  assert.deepEqual(Array.from(got.counts), want.counts, `counts ${label}`);
  want.values.forEach((v, i) => {
    if (Number.isNaN(v)) assert.ok(Number.isNaN(got.values[i]), `${label} pixel ${i}`);
    else assert.ok(Math.abs(got.values[i] - v) <= 1e-5 * Math.max(1, Math.abs(v)), `${label} pixel ${i}`);
  });
  assert.equal(got.coverage, want.counts.filter(Boolean).length / want.counts.length);
}

test('averageSlab pools unique voxels for every Laue class', () => {
  // Storage shape is (L, K, H); odd counts put a bin at the origin.
  for (const [nh, nk, nl] of [[7, 7, 5], [9, 9, 3], [5, 5, 1]]) {
    const shape = [nl, nk, nh], volume = randomVolume(shape, nh * 97 + nl);
    for (const [name] of PRESETS) {
      if (name.startsWith('m-3')) continue; // cubic groups: see the next test
      const maps = indexMaps(group(name), dims(nh, nk, nl));
      for (const fixed of [0, 1, 2]) {
        const N = [nh, nk, nl][fixed];
        for (const ids of [[0], [N >> 1], [...Array(N).keys()], [N - 1]]) {
          check(volume, shape, fixed, ids, maps, `${name} ${shape} fixed=${fixed} ids=${ids}`);
        }
      }
    }
  }
});

test('cubic groups on a cubic grid', () => {
  const shape = [5, 5, 5], volume = randomVolume(shape, 5);
  const grid = dims(5, 5, 5).map((d) => ({ ...d, edges: d.edges.map((e, i) => (i - 2.5) * 0.1) }));
  for (const name of ['m-3', 'm-3m']) {
    const maps = indexMaps(group(name), grid);
    for (const fixed of [0, 1, 2]) check(volume, shape, fixed, [1, 2, 3], maps, name);
  }
});

test('identity and even grids work without an origin bin', () => {
  const shape = [4, 6, 8], volume = randomVolume(shape, 11);
  for (const fixed of [0, 1, 2]) {
    for (const ids of [[0], [1, 2], [...Array([8, 6, 4][fixed]).keys()]]) {
      check(volume, shape, fixed, ids, [IDENTITY_MAP], `identity fixed=${fixed}`);
      // Inversion still maps bin centers onto bin centers on even grids.
      check(volume, shape, fixed, ids, indexMaps(group('-1'), dims(8, 6, 4)), `inversion fixed=${fixed}`);
    }
  }
});

test('averageSlab lays rows along the larger remaining dimension', () => {
  const shape = [2, 3, 4];
  const volume = Float32Array.from({ length: 24 }, (_, i) => i);
  const kl = averageSlab(volume, shape, 0, [1]);
  assert.deepEqual([kl.rows, kl.cols], [2, 3]);
  assert.deepEqual(Array.from(kl.values), [1, 5, 9, 13, 17, 21]);
  const hk = averageSlab(volume, shape, 2, [1]);
  assert.deepEqual([hk.rows, hk.cols], [3, 4]);
  assert.deepEqual(Array.from(hk.values), [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]);
});

test('selectBins uses bin centers within half the thickness', () => {
  const edges = Array.from({ length: 11 }, (_, i) => -0.5 + i * 0.1);
  assert.deepEqual(selectBins(edges, 0, 0.2), [4, 5]);
  assert.deepEqual(selectBins(edges, 0, 0.3), [3, 4, 5, 6]);
  assert.deepEqual(selectBins(edges, 0.05, 0.3), [4, 5, 6]);
  assert.deepEqual(selectBins(edges, 0.05, 0.1), [5]);
  assert.deepEqual(selectBins(edges, 3, 0.1), []);
});
