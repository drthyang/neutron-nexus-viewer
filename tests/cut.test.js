import assert from 'node:assert/strict';
import test from 'node:test';

import { cutAxis, cutEdges, lineCut } from '../js/cut.js';
import { IDENTITY_MAP, selectBins } from '../js/slab.js';
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

// Direct statement of the semantics: every slab voxel whose centre is within the
// band joins the bin it projects into, and a bin pools the distinct in-range
// images of its voxels under the group.
function reference(volume, shape, grid, cut, maps, variance = null) {
  const n = [shape[2], shape[1], shape[0]];
  const { fixed, x, y, a, b, dom, edges, half, geometry: { lx, ly, cos } } = cut;
  const sin = Math.sqrt(1 - cos * cos), cart = (u, v) => [lx * u + ly * cos * v, ly * sin * v];
  const A = cart(...a), B = cart(...b), L = Math.hypot(B[0] - A[0], B[1] - A[1]);
  const T = [(B[0] - A[0]) / L, (B[1] - A[1]) / L], count = edges.length - 1, step = (edges[count] - edges[0]) / count;
  const pooled = Array.from({ length: count }, () => new Set());
  for (const k0 of selectBins(grid[fixed].edges, cut.center, cut.thickness)) {
    for (let j = 0; j < n[y]; j++) {
      for (let i = 0; i < n[x]; i++) {
        const [X, Y] = cart(center(grid[x], i), center(grid[y], j)), rx = X - A[0], ry = Y - A[1];
        if (Math.abs(ry * T[0] - rx * T[1]) > half + 1e-9) continue;
        const k = Math.floor((a[dom] + (rx * T[0] + ry * T[1]) * (b[dom] - a[dom]) / L - edges[0]) / step + 1e-9);
        if (k < 0 || k >= count) continue;
        const p = [0, 0, 0];
        p[fixed] = k0; p[x] = i; p[y] = j;
        for (const { M, t } of maps) {
          const q = [0, 1, 2].map((r) => M[3 * r] * p[0] + M[3 * r + 1] * p[1] + M[3 * r + 2] * p[2] + t[r]);
          if (q.every((qi, r) => qi >= 0 && qi < n[r])) pooled[k].add(q[0] + n[0] * (q[1] + n[1] * q[2]));
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
  want.forEach(({ mean, sigma, voxels }, k) => {
    assert.equal(got.voxels[k], voxels, `${label} voxels in bin ${k}`);
    if (Number.isNaN(mean)) assert.ok(Number.isNaN(got.intensity[k]), `${label} bin ${k}`);
    else assert.ok(Math.abs(got.intensity[k] - mean) <= 1e-9 * Math.max(1, Math.abs(mean)), `${label} bin ${k}`);
    if (variance && !Number.isNaN(sigma)) assert.ok(Math.abs(got.sigma[k] - sigma) <= 1e-9 * sigma, `${label} σ in bin ${k}`);
  });
}

function makeCut(grid, fixed, a, b, { width = 0.25, thickness = 0.25, step = 0.1, geometry = { lx: 1, ly: 1, cos: 0 } } = {}) {
  const [x, y] = [0, 1, 2].filter((d) => d !== fixed);
  const dom = cutAxis(a, b);
  const [p, q] = a[dom] <= b[dom] ? [a, b] : [b, a];
  return { fixed, x, y, center: 0, thickness, a: p, b: q, dom, edges: cutEdges(p, q, dom, step), half: width / 2, geometry };
}

test('an axis-aligned cut one voxel wide reads the voxels of its row', () => {
  const shape = [5, 9, 11], grid = dims(11, 9, 5), volume = randomVolume(shape, 3);
  // Along H at K = 0.1 (row 5), L = 0 (bin 2), from H = -0.4 to 0.4: columns 1 to 9.
  const cut = makeCut(grid, 2, [-0.4, 0.1], [0.4, 0.1], { width: 0.05, thickness: 0.1 });
  const got = lineCut(volume, shape, grid, cut);
  assert.equal(got.intensity.length, 9);
  for (let k = 0; k < 9; k++) {
    const v = volume[(k + 1) + 11 * (5 + 9 * 2)];
    if (Number.isNaN(v)) assert.ok(Number.isNaN(got.intensity[k]));
    else assert.ok(Math.abs(got.intensity[k] - v) < 1e-6, `column ${k + 1}`);
    assert.equal(got.voxels[k], Number.isNaN(v) ? 0 : 1);
  }
  assert.equal(got.bins, 1);
});

test('lineCut pools unique voxels for every Laue class', () => {
  const shape = [5, 9, 9], grid = dims(9, 9, 5), volume = randomVolume(shape, 11);
  const hex = { lx: 1, ly: 1, cos: 0.5 };
  const lines = [
    [2, [-0.4, 0], [0.4, 0], {}],
    [2, [-0.3, -0.3], [0.3, 0.3], { width: 0.3 }],
    [2, [0.35, -0.2], [-0.25, 0.4], { width: 0.12, step: 0.05, geometry: hex }],
    [1, [-0.4, -0.5], [0.4, 0.5], { width: 0.2, thickness: 0.3 }],
    [0, [0, -0.5], [0.1, 0.5], { width: 0.4, step: 0.25 }],
  ];
  for (const [name] of PRESETS) {
    if (name.startsWith('m-3')) continue; // cubic groups need a cubic grid
    const maps = indexMaps(group(name), grid);
    for (const [fixed, a, b, options] of lines) {
      check(volume, shape, grid, makeCut(grid, fixed, a, b, options), maps, `${name} fixed=${fixed} ${a}→${b}`);
    }
  }
});

test('σ is that of a mean of independent voxels', () => {
  const shape = [5, 9, 9], grid = dims(9, 9, 5), volume = randomVolume(shape, 23);
  const variance = Float32Array.from(volume, (v, i) => 0.5 + (i % 7));
  const maps = indexMaps(group('6/mmm'), grid);
  check(volume, shape, grid, makeCut(grid, 2, [-0.4, 0], [0.4, 0.2], { width: 0.3 }), maps, '6/mmm', variance);
  check(volume, shape, grid, makeCut(grid, 2, [-0.4, 0], [0.4, 0.2], { width: 0.3 }), [IDENTITY_MAP], '1', variance);
});

test('the mask removes voxels, and inverted keeps only those', () => {
  const shape = [3, 5, 7], grid = dims(7, 5, 3), volume = Float32Array.from({ length: 105 }, (_, i) => i);
  const mask = Uint8Array.from({ length: 105 }, (_, i) => (i % 2));
  const cut = makeCut(grid, 2, [-0.3, 0], [0.3, 0], { width: 0.05, thickness: 0.1 });
  const kept = lineCut(volume, shape, grid, cut, [IDENTITY_MAP], mask, false);
  const removed = lineCut(volume, shape, grid, cut, [IDENTITY_MAP], mask, true);
  // Row K = 0 (2) of the middle L bin (1): flat indices 7 * (2 + 5) + i.
  for (let k = 0; k < 7; k++) {
    const f = 49 + k, odd = f % 2 === 1;
    assert.equal(Number.isNaN(kept.intensity[k]), odd);
    assert.equal(Number.isNaN(removed.intensity[k]), !odd);
  }
});

test('cutEdges centres bins on the start and steps to the end', () => {
  const edges = cutEdges([-0.3, 0], [0.3, 0.1], 0, 0.1);
  assert.equal(edges.length, 8);
  assert.ok(Math.abs(edges[0] + 0.35) < 1e-12 && Math.abs(edges[7] - 0.35) < 1e-12);
  assert.equal(cutAxis([0, 0], [0.1, 0.5]), 1);
  assert.equal(cutAxis([0, 0], [0.5, 0.5]), 0);
  assert.throws(() => cutEdges([0, 0], [1, 0], 0, 1e-6), /too fine/);
});
