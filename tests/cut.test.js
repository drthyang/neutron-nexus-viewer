import assert from 'node:assert/strict';
import test from 'node:test';

import { cutAxis, cutBand, cutEdges, lineCut } from '../js/cut.js';
import { IDENTITY_MAP } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps, PRESETS } from '../js/symmetry.js';

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
const center = (d, i) => d.edges[0] + (i + 0.5) * (d.edges[1] - d.edges[0]);
const UNIT = [1, 0, 0, 0, 1, 0, 0, 0, 1];
// A hexagonal reciprocal basis (a* and b* at 60°, c* shorter), row-major display -> Cartesian.
const HEX = [1, 0.5, 0, 0, Math.sqrt(3) / 2, 0, 0, 0, 0.6];
const cart = (T, x) => [0, 1, 2].map((i) => T[3 * i] * x[0] + T[3 * i + 1] * x[1] + T[3 * i + 2] * x[2]);

// Direct statement of the semantics: every voxel whose centre lies within the
// radius of the line (in 3-D) joins the bin it projects into, and a bin pools
// the distinct in-range images of its voxels under the group.
function reference(volume, shape, grid, cut, maps, variance = null) {
  const n = [shape[2], shape[1], shape[0]];
  const { a, b, dom, edges, radius, T } = cut;
  const A = cart(T, a), B = cart(T, b), L = Math.hypot(...B.map((v, i) => v - A[i])), dir = B.map((v, i) => (v - A[i]) / L);
  const count = edges.length - 1, step = (edges[count] - edges[0]) / count;
  const pooled = Array.from({ length: count }, () => new Set());
  for (let k2 = 0; k2 < n[2]; k2++) {
    for (let k1 = 0; k1 < n[1]; k1++) {
      for (let k0 = 0; k0 < n[0]; k0++) {
        const p = [k0, k1, k2], X = cart(T, p.map((i, d) => center(grid[d], i))), r = X.map((v, i) => v - A[i]);
        const t = r.reduce((s, v, i) => s + v * dir[i], 0);
        if (Math.hypot(...r.map((v, i) => v - t * dir[i])) > radius * (1 + 1e-9)) continue;
        const k = Math.floor((a[dom] + t * (b[dom] - a[dom]) / L - edges[0]) / step + 1e-9);
        if (k < 0 || k >= count) continue;
        for (const { M, t: shift } of maps) {
          const q = [0, 1, 2].map((i) => M[3 * i] * p[0] + M[3 * i + 1] * p[1] + M[3 * i + 2] * p[2] + shift[i]);
          if (q.every((qi, i) => qi >= 0 && qi < n[i])) pooled[k].add(q[0] + n[0] * (q[1] + n[1] * q[2]));
        }
      }
    }
  }
  return pooled.map((set) => {
    const finite = [...set].filter((f) => Number.isFinite(volume[f]));
    const mean = finite.length ? finite.reduce((s, f) => s + volume[f], 0) / finite.length : NaN;
    const sigma = finite.length && variance ? Math.sqrt(finite.reduce((s, f) => s + variance[f], 0)) / finite.length : NaN;
    return { mean, sigma, voxels: finite.length };
  });
}

function check(volume, shape, grid, cut, maps, label, variance = null) {
  const got = lineCut(volume, shape, grid, cut, maps, null, false, variance);
  const want = reference(volume, shape, grid, cut, maps, variance);
  assert.ok(want.some((w) => w.voxels > 0), `${label}: the cut should meet data`);
  want.forEach(({ mean, sigma, voxels }, k) => {
    assert.equal(got.voxels[k], voxels, `${label} voxels in bin ${k}`);
    if (Number.isNaN(mean)) assert.ok(Number.isNaN(got.intensity[k]), `${label} bin ${k}`);
    else assert.ok(Math.abs(got.intensity[k] - mean) <= 1e-9 * Math.max(1, Math.abs(mean)), `${label} bin ${k}`);
    if (variance && !Number.isNaN(sigma)) assert.ok(Math.abs(got.sigma[k] - sigma) <= 1e-9 * sigma, `${label} σ in bin ${k}`);
  });
}

function makeCut(a, b, { width = 0.25, step = 0.1, T = UNIT } = {}) {
  const dom = cutAxis(a, b), [p, q] = a[dom] <= b[dom] ? [a, b] : [b, a];
  return { a: p, b: q, dom, edges: cutEdges(p, q, dom, step), radius: width / 2, T };
}

test('a thin rod along an axis reads the voxels on its line', () => {
  const shape = [5, 9, 11], grid = dims(11, 9, 5), volume = randomVolume(shape, 3);
  // Along H at K = 0.1 (row 5), L = 0 (bin 2), from H = -0.4 to 0.4: columns 1 to 9.
  const got = lineCut(volume, shape, grid, makeCut([-0.4, 0.1, 0], [0.4, 0.1, 0], { width: 0.05 }));
  assert.equal(got.intensity.length, 9);
  for (let k = 0; k < 9; k++) {
    const v = volume[(k + 1) + 11 * (5 + 9 * 2)];
    if (Number.isNaN(v)) assert.ok(Number.isNaN(got.intensity[k]));
    else assert.ok(Math.abs(got.intensity[k] - v) < 1e-6, `column ${k + 1}`);
    assert.equal(got.voxels[k], Number.isNaN(v) ? 0 : 1);
  }
});

test('W widens the rod across the line in every direction', () => {
  // Along H through voxel centres, a rod of diameter 2.2 bins takes the row
  // and its four neighbours across (in K and in L): five voxels per point.
  const shape = [5, 5, 9], grid = dims(9, 5, 5).map((d) => ({ ...d, edges: d.edges.map((e, i) => (i - (d.edges.length - 1) / 2) * 0.1) }));
  const volume = new Float32Array(9 * 5 * 5).fill(1);
  const got = lineCut(volume, shape, grid, makeCut([-0.3, 0, 0], [0.3, 0, 0], { width: 0.22 }));
  assert.deepEqual(Array.from(got.voxels), [5, 5, 5, 5, 5, 5, 5]);
});

test('lineCut pools unique voxels for every Laue class', () => {
  const shape = [5, 9, 9], grid = dims(9, 9, 5), volume = randomVolume(shape, 11);
  const lines = [
    [[-0.4, 0, 0], [0.4, 0, 0], {}],
    [[-0.3, -0.3, 0], [0.3, 0.3, 0], { width: 0.3 }],
    [[0.35, -0.2, 0], [-0.25, 0.4, 0], { width: 0.12, step: 0.05, T: HEX }],
    [[-0.4, 0, -0.5], [0.4, 0, 0.5], { width: 0.2, T: HEX }],
    [[0, -0.4, -0.5], [0, 0.1, 0.5], { width: 0.4, step: 0.25 }],
  ];
  for (const [name] of PRESETS) {
    if (name.startsWith('m-3')) continue; // cubic groups need a cubic grid
    const maps = indexMaps(group(name), grid);
    for (const [a, b, options] of lines) check(volume, shape, grid, makeCut(a, b, options), maps, `${name} ${a}→${b}`);
  }
});

test('σ is that of a mean of independent voxels', () => {
  const shape = [5, 9, 9], grid = dims(9, 9, 5), volume = randomVolume(shape, 23);
  const variance = Float32Array.from(volume, (v, i) => 0.5 + (i % 7));
  const cut = makeCut([-0.4, 0, 0], [0.4, 0.2, 0], { width: 0.3, T: HEX });
  check(volume, shape, grid, cut, indexMaps(group('6/mmm'), grid), '6/mmm', variance);
  check(volume, shape, grid, cut, [IDENTITY_MAP], '1', variance);
});

test('the mask removes voxels, and inverted keeps only those', () => {
  const shape = [3, 5, 7], grid = dims(7, 5, 3), volume = Float32Array.from({ length: 105 }, (_, i) => i);
  const mask = Uint8Array.from({ length: 105 }, (_, i) => (i % 2));
  const cut = makeCut([-0.3, 0, 0], [0.3, 0, 0], { width: 0.05 });
  const kept = lineCut(volume, shape, grid, cut, [IDENTITY_MAP], mask, false);
  const removed = lineCut(volume, shape, grid, cut, [IDENTITY_MAP], mask, true);
  // Row K = 0 (2) of the middle L bin (1): flat indices 7 * (2 + 5) + i.
  for (let k = 0; k < 7; k++) {
    const f = 49 + k, odd = f % 2 === 1;
    assert.equal(Number.isNaN(kept.intensity[k]), odd);
    assert.equal(Number.isNaN(removed.intensity[k]), !odd);
  }
});

test('the band is the rod across its slice, W wide', () => {
  const cut = makeCut([-0.3, 0.1, 0], [0.3, 0.1, 0], { width: 0.2, T: HEX });
  const band = cutBand(cut, 2);
  // In the plane L = 0, from the first bin edge to the last, 0.1 to either side in length.
  for (const corner of band) assert.ok(Math.abs(corner[2]) < 1e-12);
  const across = (c) => { const X = cart(HEX, c), A = cart(HEX, cut.a); return Math.abs(X[1] - A[1]); };
  for (const corner of band) assert.ok(Math.abs(across(corner) - 0.1) < 1e-12);
  // The first two corners straddle the line at its first bin edge.
  const start = [0, 1, 2].map((d) => (band[0][d] + band[3][d]) / 2);
  assert.ok(Math.abs(start[0] - cut.edges[0]) < 1e-12 && Math.abs(start[1] - 0.1) < 1e-12);
});

test('cutEdges centres bins on the start and steps to the end', () => {
  const edges = cutEdges([-0.3, 0], [0.3, 0.1], 0, 0.1);
  assert.equal(edges.length, 8);
  assert.ok(Math.abs(edges[0] + 0.35) < 1e-12 && Math.abs(edges[7] - 0.35) < 1e-12);
  assert.equal(cutAxis([0, 0], [0.1, 0.5]), 1);
  assert.equal(cutAxis([0, 0], [0.5, 0.5]), 0);
  assert.equal(cutAxis([0, 0, 0], [0.1, 0.2, 0.3]), 2);
  assert.throws(() => cutEdges([0, 0], [1, 0], 0, 1e-6), /too fine/);
});
