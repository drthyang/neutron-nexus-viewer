import assert from 'node:assert/strict';
import { unlinkSync } from 'node:fs';
import { tmpdir } from 'node:os';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { describeFile, loadVolume } from '../js/nexus.js';
import { formatIsawUB, parseIsawUB } from '../js/reindex.js';
import { bMatrix, cellFromUB, inv3, levenbergMarquardt, mul, mulv, rot } from '../js/rigaku-geometry.js';
import { writeMantidMD } from '../js/savemd.js';
import {
  UB_SYSTEMS, findPeaks, fitsSystem, formatSU, indexShift, peaksCSV, refineUB, startParams, suggestSystem, voxelQ,
} from '../js/ub-refine.js';

await h5wasm.ready;
const example = (name) => new URL(`../examples/${name}`, import.meta.url).pathname;
const HEX = [4.2, 4.2, 6.8, 90, 90, 120];

/** Display axes [H,0,0], [0,K,0], [0,0,L] with the same bin edges. */
const hklDims = (edges) => ['[H,0,0]', '[0,K,0]', '[0,0,L]'].map((longName, d) => ({
  longName, label: 'HKL'[d], frame: 'HKL', units: 'r.l.u.', edges, basis: { vec: [0, 1, 2].map((i) => (i === d ? 1 : 0)), letter: 'HKL'[d] },
}));

/**
 * A hexagonal crystal (UB `ubTrue`) binned with a wrong UB `ubOld`, as a reduction with an
 * error in its UB would bin it: Gaussian Bragg peaks of σ = `sigma` Å⁻¹ at Q = ubTrue·h, on a
 * background of 1, with Gaussian noise of `noise` × √I (a fixed seed).
 */
function binnedCrystal(ubTrue, ubOld, { step = 0.05, sigma = 0.012, noise = 0 } = {}) {
  const half = 4 + step / 2, nb = Math.round(2 * half / step);
  const edges = Array.from({ length: nb + 1 }, (_, k) => -half + k * step), c = edges.slice(0, -1).map((e) => e + step / 2);
  const inv = inv3(ubTrue), volume = new Float32Array(nb ** 3);
  let seed = 1;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let i0 = 0, v = 0; i0 < nb; i0++) {
    for (let i1 = 0; i1 < nb; i1++) {
      for (let i2 = 0; i2 < nb; i2++, v++) {
        const q = mulv(ubOld, [c[i2], c[i1], c[i0]]), n = mulv(inv, q).map(Math.round);
        let value = 1;
        if (n.some(Boolean)) {
          const Q = mulv(ubTrue, n), d2 = (q[0] - Q[0]) ** 2 + (q[1] - Q[1]) ** 2 + (q[2] - Q[2]) ** 2;
          value += 1000 * Math.exp(-(Q[0] ** 2 + Q[1] ** 2 + Q[2] ** 2)) * Math.exp(-d2 / (2 * sigma * sigma));
        }
        volume[v] = noise ? value + noise * Math.sqrt(value) * Math.sqrt(-2 * Math.log(rnd())) * Math.cos(2 * Math.PI * rnd()) : value;
      }
    }
  }
  return { volume, shape: [nb, nb, nb], dims: hklDims(edges) };
}

const U = rot([1, 2, 3], 25);
const ubTrue = mul(U, bMatrix(...HEX));
// 0.2-0.4 % off in every length, 0.2-0.4° in the angles, and turned by 0.4°.
const ubOld = mul(mul(rot([0.3, -1, 0.5], 0.4), U), bMatrix(4.21, 4.19, 6.83, 90.3, 89.8, 120.4));

for (const noise of [0, 0.5]) {
  test(`refinement recovers the UB a volume was binned with a wrong one${noise ? ', with noise' : ''}`, () => {
    const { volume, shape, dims } = binnedCrystal(ubTrue, ubOld, { noise });
    const { peaks, stats } = findPeaks(volume, shape, dims);
    assert.equal(stats.radius, 5, 'box from half widths of about 1 and 2 voxels');
    assert.ok(peaks.length > 300, `${peaks.length} peaks`);
    const scale = Math.cbrt(voxelQ(dims, ubOld).reduce((p, x) => p * x, 1));
    const hex = refineUB(peaks, ubOld, 'hexagonal', { scale });
    assert.ok(hex.settled && hex.used === peaks.length, `${hex.used} of ${peaks.length} used`);
    // The cell to 12 ppm, against an error of 0.2-0.4 % in the UB the volume was binned with.
    assert.ok(Math.abs(hex.cell[0] - 4.2) < 5e-5 && Math.abs(hex.cell[2] - 6.8) < 1e-4, `cell ${hex.cell}`);
    assert.deepEqual(hex.cell.slice(3), [90, 90, 120]);
    assert.equal(hex.cell[0], hex.cell[1]);
    assert.ok(hex.ub.every((x, i) => Math.abs(x - ubTrue[i]) < 1e-5), 'UB elements');
    assert.ok(hex.su[0] > 0 && hex.su[2] > 0 && hex.su.slice(3).every((x) => x === 0) && hex.volumeSU > 0);
    const tri = refineUB(peaks, ubOld, 'triclinic', { scale });
    assert.ok(tri.cell.slice(0, 3).every((x, i) => Math.abs(x - HEX[i]) < 1e-4) && tri.cell.slice(3).every((x, i) => Math.abs(x - HEX[3 + i]) < 2e-3), `triclinic ${tri.cell}`);
    // Reindexing with the refined UB moves the corner indices (±4) by the old UB's error.
    const expected = indexShift(mul(inv3(ubTrue), ubOld), dims);
    assert.ok(Math.abs(indexShift(mul(inv3(hex.ub), ubOld), dims) - expected) < 1e-3);
    const csv = peaksCSV(peaks, hex).trim().split('\n');
    assert.equal(csv.length, peaks.length + 1);
    assert.equal(csv[0], 'H,K,L,h_volume,k_volume,l_volume,h_refined,k_refined,l_refined,d_A,residual_invA_no2pi,max,intensity,used');
  });
}

test('peaks off the lattice are left out by their residual', () => {
  const { volume, shape, dims } = binnedCrystal(ubTrue, ubTrue);
  const { peaks } = findPeaks(volume, shape, dims);
  // Displace five peaks by 0.08 r.l.u. along l: within the index tolerance, far off the fit.
  const moved = peaks.map((p, k) => (k % 50 === 7 ? { ...p, h: [p.h[0], p.h[1], p.h[2] + 0.08] } : p));
  const scale = Math.cbrt(voxelQ(dims, ubTrue).reduce((p, x) => p * x, 1));
  const r = refineUB(moved, ubTrue, 'hexagonal', { scale });
  const off = moved.filter((_, k) => k % 50 === 7).length;
  assert.equal(r.within, moved.length);
  assert.equal(r.used, moved.length - off);
  assert.ok(Math.abs(r.cell[0] - 4.2) < 1e-5 && Math.abs(r.cell[2] - 6.8) < 1e-5, `cell ${r.cell}`);
});

test('the example: the refined cell is the one it was generated with', () => {
  const file = new h5wasm.File(example('demo_10K.nxs'), 'r');
  let info, volume;
  try {
    info = describeFile(file);
    ({ volume } = loadVolume(file, info));
  } finally {
    file.close();
  }
  const { peaks } = findPeaks(volume, info.shape, info.dims);
  const scale = Math.cbrt(voxelQ(info.dims, info.lattice.ub).reduce((p, x) => p * x, 1));
  assert.equal(suggestSystem(cellFromUB(info.lattice.ub)).system, 'hexagonal');
  const r = refineUB(peaks, info.lattice.ub, 'hexagonal', { scale });
  // a = 4.2, c = 6.8 and the file's UB exact: superlattice peaks (at the M points) are not indexed.
  assert.ok(r.used > 200 && r.used < peaks.length, `${r.used} of ${peaks.length}`);
  assert.ok(Math.abs(r.cell[0] - 4.2) < 1e-4 && Math.abs(r.cell[2] - 6.8) < 2e-4, `cell ${r.cell}`);
  assert.ok(indexShift(mul(inv3(r.ub), info.lattice.ub), info.dims) < 2e-3);
});

test('crystal systems: constraints, suggestions and refusals', () => {
  assert.deepEqual(UB_SYSTEMS.rhombohedral.cell(startParams('rhombohedral', [5, 5.1, 4.9, 80, 81, 79])), [5, 5, 5, 80, 80, 80]);
  assert.deepEqual(UB_SYSTEMS['monoclinic-c'].cell(startParams('monoclinic-c', [5, 6, 7, 90, 90, 100])), [5, 6, 7, 90, 90, 100]);
  assert.ok(fitsSystem('hexagonal', [4.2, 4.21, 6.8, 90.2, 89.9, 120.3]));
  assert.ok(!fitsSystem('hexagonal', [4.2, 4.2, 6.8, 90, 90, 60]));
  assert.deepEqual(suggestSystem(HEX, '6/mmm'), { system: 'hexagonal', from: 'laue' });
  assert.deepEqual(suggestSystem(HEX, '-3m1'), { system: 'hexagonal', from: 'laue' });
  assert.deepEqual(suggestSystem(HEX, '4/mmm'), { system: 'hexagonal', from: 'metric' }); // the Laue class does not fit the cell
  assert.deepEqual(suggestSystem([5, 5, 5, 90, 90, 90], 'm-3m'), { system: 'cubic', from: 'laue' });
  assert.deepEqual(suggestSystem([5, 5, 5, 90, 90, 90], 'mmm'), { system: 'orthorhombic', from: 'laue' });
  assert.equal(suggestSystem([5, 6, 7, 90, 90, 90]).system, 'orthorhombic');
  assert.equal(suggestSystem([5, 6, 7, 90, 95, 90]).system, 'monoclinic');
  assert.equal(suggestSystem([5, 6, 7, 90, 90, 95], '2/m (c unique)').system, 'monoclinic-c');
  assert.equal(suggestSystem([5, 6, 7, 91.5, 95, 90]).system, 'triclinic');
  assert.throws(() => refineUB([], ubTrue, 'cubic'), /does not have the metric of the cubic system/);
  assert.throws(() => refineUB([], mul(ubTrue, [1, 0, 0, 0, 1, 0, 0, 0, -1]), 'triclinic'), /left-handed/);
});

test('standard uncertainties are written as crystallographers write them', () => {
  assert.equal(formatSU(4.200031, 2.1e-5), '4.20003(2)');
  assert.equal(formatSU(6.80005, 6.6e-5), '6.80005(7)');
  assert.equal(formatSU(10.02738, 0.0015), '10.0274(15)');
  assert.equal(formatSU(10.02738, 0.0025), '10.027(3)');
  assert.equal(formatSU(554.3, 12), '554(12)');
  assert.equal(formatSU(1234.5, 25), '1235(25)');
  assert.equal(formatSU(90, 0), '90.0000');
  assert.equal(formatSU(90, NaN, 2), '90.00');
});

test('the least squares returns the covariance of the parameters', () => {
  const xs = Array.from({ length: 20 }, (_, i) => i / 4), ys = xs.map((x, i) => 2 * x + 1 + 0.01 * Math.sin(7 * i));
  const fit = levenbergMarquardt((p) => Float64Array.from(xs, (x, i) => p[0] * x + p[1] - ys[i]), [1, 0], { steps: [1e-6, 1e-6], fScale: 1e6 });
  const C = fit.covariance;
  assert.ok(Math.abs(C[1] - C[2]) < 1e-15 && C[1] < 0, 'symmetric, slope and intercept anticorrelated');
  assert.ok(Math.abs(Math.sqrt(C[0]) - fit.sigma[0]) < 1e-15 && Math.abs(Math.sqrt(C[3]) - fit.sigma[1]) < 1e-15);
});

test('an ISAW UB file round-trips the UB, the cell and its uncertainties', () => {
  const errors = { a: 2e-5, b: 2e-5, c: 7e-5, alpha: 0, beta: 0, gamma: 0, volume: 1.3e-3 };
  const text = formatIsawUB(ubTrue, { errors, note: 'UB refined in a test.' });
  const parsed = parseIsawUB(text);
  assert.ok(parsed.ub.every((x, i) => Math.abs(x - ubTrue[i]) < 1e-10));
  assert.ok(['a', 'b', 'c', 'alpha', 'beta', 'gamma'].every((k, i) => Math.abs(parsed.cell[k] - HEX[i]) < 1e-6));
  assert.deepEqual(parsed.errors, errors);
  assert.match(text, /\n\nUB refined in a test\.\n/);
  assert.equal(parseIsawUB(formatIsawUB(ubTrue)).errors.a, 0);
});

test('the Mantid writer stores the cell\'s uncertainties', () => {
  const path = `${tmpdir()}/ub-refine-errors-${process.pid}.nxs`, edges = [[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]];
  writeMantidMD(h5wasm, path, {
    shape: [2, 2, 2], edges, signal: new Float64Array(8), errors2: new Float64Array(8), events: new Float64Array(8),
    ub: ubTrue, cell: HEX, cellErrors: [2e-5, 2e-5, 7e-5, 0, 0, 0],
  });
  const file = new h5wasm.File(path, 'r');
  try {
    const lattice = file.get('MDHistoWorkspace/experiment0/sample/oriented_lattice');
    assert.deepEqual(['a', 'c', 'gamma'].map((k) => lattice.get(`unit_cell_${k}_error`).value[0]), [2e-5, 7e-5, 0]);
  } finally {
    file.close();
    unlinkSync(path);
  }
});
