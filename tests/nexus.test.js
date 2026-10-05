import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { cartesianBasis, cellFromUB, describeFile, loadVariance, loadVolume, nominalCell, parseBasis, planeGeometry, reciprocalMetric } from '../js/nexus.js';
import { powderPlan } from '../js/powder.js';
import { averageSlab, IDENTITY_MAP, selectBins } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps } from '../js/symmetry.js';

const fixture = (name) => new URL(`./fixtures/${name}`, import.meta.url).pathname;
const expected = JSON.parse(readFileSync(fixture('expected.json')));
await h5wasm.ready;

function open(name, options) {
  const file = new h5wasm.File(fixture(name), 'r');
  try {
    const info = describeFile(file, options);
    return { info, ...loadVolume(file, info), variance: loadVariance(file, info) };
  } finally {
    file.close();
  }
}

function checkCases({ info, volume }, cases) {
  for (const c of cases) {
    const label = `fixed=${c.fixed} center=${c.center} thickness=${c.thickness} inversion=${c.inversion}`;
    const ids = selectBins(info.dims[c.fixed].edges, c.center, c.thickness);
    assert.deepEqual(ids, c.ids, label);
    const maps = c.inversion ? indexMaps(closeGroup(parseOps('-h,-k,-l')), info.dims) : [IDENTITY_MAP];
    const got = averageSlab(volume, info.shape, c.fixed, ids, maps);
    assert.deepEqual([got.rows, got.cols], c.shape, label);
    assert.deepEqual(Array.from(got.counts), c.counts, label);
    c.values.forEach((v, i) => {
      if (v === null) assert.ok(Number.isNaN(got.values[i]), `${label} pixel ${i}`);
      else assert.ok(Math.abs(got.values[i] - v) <= 1e-5 * Math.max(1, Math.abs(v)), `${label} pixel ${i}: ${got.values[i]} vs ${v}`);
    });
  }
}

test('Mantid MDHistoWorkspace fixture matches the Python reference', () => {
  const loaded = open('mdhisto_small.nxs');
  const { info } = loaded;
  assert.equal(info.signal, '/MDHistoWorkspace/data/signal');
  assert.deepEqual(info.mask, { path: '/MDHistoWorkspace/data/mask', valid: false });
  assert.deepEqual(info.shape, [7, 9, 11]);
  assert.deepEqual(info.dims.map((d) => d.label), ['H', 'K', 'L']);
  assert.deepEqual(info.dims.map((d) => d.edges.length), [12, 10, 8]);
  assert.equal(info.dims[0].units, 'r.l.u.');
  // The cell comes from the UB matrix and agrees with the stored unit_cell values.
  assert.equal(info.lattice.source, 'UB');
  for (const [k, v] of Object.entries({ a: 8.03, b: 8.02, c: 10.03, alpha: 90.1, beta: 90.2, gamma: 119.9 })) {
    assert.ok(Math.abs(info.lattice[k] - v) < 1e-9, k);
  }
  checkCases(loaded, expected.mdhisto_small);
  // errors_squared holds |signal| in this fixture; non-finite entries are unknown.
  assert.deepEqual(info.errors, { path: '/MDHistoWorkspace/data/errors_squared', squared: true });
  let known = 0;
  loaded.volume.forEach((v, i) => {
    if (Number.isNaN(v)) return;
    assert.equal(loaded.variance[i], Math.abs(v), `voxel ${i}`);
    known++;
  });
  assert.ok(known > 500);
  assert.ok(loaded.variance.some(Number.isNaN));
});

test('plain NXdata fixture: size-1 axis dropped, centers become edges', () => {
  const loaded = open('nxdata_small.h5');
  const { info } = loaded;
  assert.equal(info.signal, '/entry/data/intensity');
  assert.equal(info.mask, null);
  assert.equal(info.errors, null);
  assert.equal(loaded.variance, null);
  assert.deepEqual(info.shape, [6, 8, 10]);
  assert.deepEqual(info.dims.map((d) => d.name), ['qx', 'qy', 'qz']);
  assert.ok(Math.abs(info.dims[2].edges[0] - 0.25) < 1e-12);
  assert.equal(info.lattice, null);
  const geometry = planeGeometry(info.dims, info.lattice, 0, 1);
  assert.deepEqual(geometry, { lx: 1, ly: 1, cos: 0, equal: true, lattice: false });
  // Axes in Å⁻¹ give |Q| without a cell.
  assert.equal(powderPlan(info.dims, info.lattice).frame, 'Q');
  checkCases(loaded, expected.nxdata_small);
});

test('NEBULA3D volume: (H, K, L) storage, mask 1 = valid, sigma, UB with 2π', () => {
  const loaded = open('nebula3d_small.h5');
  const { info } = loaded;
  assert.equal(info.signal, '/entry/data');
  assert.deepEqual(info.mask, { path: '/entry/mask', valid: true });
  assert.deepEqual(info.errors, { path: '/entry/sigma', squared: false });
  // Displayed fastest storage axis first, as for Mantid's D0, D1, D2.
  assert.deepEqual(info.shape, [5, 7, 9]);
  assert.deepEqual(info.dims.map((d) => d.label), ['L', 'K', 'H']);
  assert.ok(info.dims.every((d) => d.units === 'r.l.u.' && d.frame === 'HKL'));
  assert.ok(Math.abs(info.dims[2].edges[0] + 0.5) < 1e-12);
  assert.equal(info.signed, false);
  // ub_matrix carries 2π; the cell is that of the UB without it.
  assert.equal(info.lattice.source, 'UB');
  for (const [k, v] of Object.entries({ a: 5.91, b: 10.42, c: 24.79, alpha: 89.55, beta: 90.61, gamma: 90.63 })) {
    assert.ok(Math.abs(info.lattice[k] - v) < 1e-9, k);
  }
  checkCases(loaded, expected.nebula3d_small);
  // σ is squared when read.
  let known = 0;
  loaded.volume.forEach((v, i) => {
    if (Number.isNaN(v)) return;
    assert.ok(loaded.variance[i] > 0.0001 && loaded.variance[i] < 1.5, `voxel ${i}`);
    known++;
  });
  assert.equal(known, loaded.stats.valid);
});

test('NEBULA3D 3D-ΔPDF: real-space axes along a, b, c, cell from lat_*, signed', () => {
  const loaded = open('nebula3d_dpdf_small.h5');
  const { info } = loaded;
  assert.equal(info.signal, '/data');
  assert.equal(info.mask, null);
  assert.equal(info.errors, null);
  assert.equal(loaded.variance, null);
  assert.equal(info.signed, true);
  assert.deepEqual(info.dims.map((d) => d.label), ['z', 'y', 'x']);
  assert.deepEqual(info.dims.map((d) => [d.frame, d.units, d.axis, d.length]), [['direct', 'Å', 2, 6], ['direct', 'Å', 1, 4], ['direct', 'Å', 0, 4]]);
  assert.deepEqual([info.lattice.a, info.lattice.gamma, info.lattice.source], [4, 120, 'lat_* attributes']);
  checkCases(loaded, expected.nebula3d_dpdf_small);
  // The x-y section is drawn at γ = 120°, the others at 90°.
  assert.ok(Math.abs(planeGeometry(info.dims, info.lattice, 2, 1).cos + 0.5) < 1e-12);
  assert.ok(Math.abs(planeGeometry(info.dims, info.lattice, 0, 1).cos) < 1e-12);
  const { T, lattice } = cartesianBasis(info.dims, info.lattice);
  assert.equal(lattice, false);
  const column = (d) => [T[d], T[3 + d], T[6 + d]];
  assert.ok(Math.abs(column(2).reduce((s, v, i) => s + v * column(1)[i], 0) + 0.5) < 1e-12);
  // The 6-fold axis, h+k,-h,l on reciprocal coordinates, is x-y,x,z on these:
  // it sends (x, y, z) = (0.5, 0, 0) to (0.5, 0.5, 0).
  const maps = indexMaps(closeGroup(parseOps('h+k,-h,l')), info.dims);
  assert.equal(maps.length, 6);
  const index = (x, y, z) => [z / 0.6 + 3, y / 0.5 + 4, x / 0.5 + 4].map(Math.round);
  const apply = ({ M, t }, i) => [0, 1, 2].map((d) => M[3 * d] * i[0] + M[3 * d + 1] * i[1] + M[3 * d + 2] * i[2] + t[d]);
  const sixfold = indexMaps(parseOps('h+k,-h,l'), info.dims)[0];
  assert.deepEqual(apply(sixfold, index(0.5, 0, 0)), index(0.5, 0.5, 0));
  assert.ok(maps.every((m) => String(apply(m, index(0, 0, 0))) === String(index(0, 0, 0))));
});

test('3D-ΔPDF in the Mantid layout (NEBULA3D): x, y, z in Å along a, b, c, zero errors count as none', () => {
  const loaded = open('mdhisto_dpdf_small.nxs');
  const { info } = loaded;
  assert.deepEqual(info.dims.map((d) => d.label), ['z', 'y', 'x']);
  assert.deepEqual(info.dims.map((d) => [d.frame, d.units, d.axis]), [['direct', 'Å', 2], ['direct', 'Å', 1], ['direct', 'Å', 0]]);
  assert.equal(info.signed, true);
  assert.equal(info.lattice.source, 'UB');
  assert.ok(Math.abs(info.lattice.gamma - 120) < 1e-9 && Math.abs(info.dims[0].length - 6) < 1e-9);
  assert.ok(Math.abs(planeGeometry(info.dims, info.lattice, 2, 1).cos + 0.5) < 1e-9);
  // errors_squared is there (Mantid needs it) but all zero: no uncertainties.
  assert.ok(info.errors);
  assert.equal(loaded.variance, null);
  checkCases(loaded, expected.mdhisto_dpdf_small);
});

test('a dataset opened with another\'s axis order is read transposed to it', () => {
  const plain = open('mdhisto_small.nxs');
  const aligned = open('mdhisto_small.nxs', { order: ['L', 'K', 'H'] });
  assert.deepEqual(aligned.info.dims.map((d) => d.label), ['L', 'K', 'H']);
  assert.deepEqual(aligned.info.shape, [11, 9, 7]);
  assert.deepEqual(aligned.stats, plain.stats);
  // Storage is (L, K, H); the aligned volume is held as (H, K, L).
  for (let l = 0; l < 7; l++) {
    for (let k = 0; k < 9; k++) {
      for (let h = 0; h < 11; h++) {
        const a = aligned.volume[(h * 9 + k) * 7 + l], p = plain.volume[(l * 9 + k) * 11 + h];
        assert.ok(Object.is(a, p) || a === p, `${h} ${k} ${l}`);
        assert.ok(Object.is(aligned.variance[(h * 9 + k) * 7 + l], plain.variance[(l * 9 + k) * 11 + h]));
      }
    }
  }
  // Labels that are not a reordering of the file's axes leave it as it is.
  assert.deepEqual(open('mdhisto_small.nxs', { order: ['x', 'y', 'z'] }).info.dims.map((d) => d.label), ['H', 'K', 'L']);
});

test('parseBasis reads Mantid axis names', () => {
  assert.deepEqual(parseBasis('[H,0,0]'), { vec: [1, 0, 0], letter: 'H' });
  assert.deepEqual(parseBasis('[-K,K,0]'), { vec: [-1, 1, 0], letter: 'K' });
  assert.deepEqual(parseBasis('[0.5H, 0.5H, 0]'), { vec: [0.5, 0.5, 0], letter: 'H' });
  assert.equal(parseBasis('DeltaE'), null);
  assert.equal(parseBasis('[H,K,0]'), null);
});

test('hexagonal reciprocal geometry', () => {
  const cell = { a: 8, b: 8, c: 10, alpha: 90, beta: 90, gamma: 120 };
  const G = reciprocalMetric(cell);
  const astar = 1 / (8 * Math.sin(Math.PI / 3));
  assert.ok(Math.abs(Math.sqrt(G[0][0]) - astar) < 1e-12);
  assert.ok(Math.abs(Math.sqrt(G[2][2]) - 0.1) < 1e-12);
  const dims = ['[H,0,0]', '[0,K,0]', '[0,0,L]'].map((n) => ({ frame: 'HKL', units: 'r.l.u.', basis: parseBasis(n) }));
  const hk = planeGeometry(dims, cell, 0, 1);
  assert.ok(Math.abs(hk.cos - 0.5) < 1e-12);
  assert.ok(Math.abs(planeGeometry(dims, cell, 0, 2).cos) < 1e-12);
  const snapped = nominalCell({ ...cell, alpha: 90.4, gamma: 119.2, beta: 103.5 });
  assert.deepEqual([snapped.alpha, snapped.beta, snapped.gamma], [90, 103.5, 120]);
});

test('cellFromUB recovers the cell from a Busing-Levy B matrix', () => {
  // Hexagonal a = 4, c = 7: a* = 1/(a sin 60), c* = 1/c, gamma* = 60 degrees.
  const astar = 1 / (4 * Math.sin(Math.PI / 3));
  const B = [astar, astar * 0.5, 0, 0, astar * Math.sin(Math.PI / 3), 0, 0, 0, 1 / 7];
  const cell = cellFromUB(B);
  for (const [k, v] of Object.entries({ a: 4, b: 4, c: 7, alpha: 90, beta: 90, gamma: 120 })) {
    assert.ok(Math.abs(cell[k] - v) < 1e-9, k);
  }
  assert.equal(cellFromUB([1, 2, 3, 2, 4, 6, 0, 0, 1]), null);
});
