import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { cellFromUB, describeFile, loadVariance, loadVolume, nominalCell, parseBasis, planeGeometry, reciprocalMetric } from '../js/nexus.js';
import { powderPlan } from '../js/powder.js';
import { averageSlab, IDENTITY_MAP, selectBins } from '../js/slab.js';
import { closeGroup, indexMaps, parseOps } from '../js/symmetry.js';

const fixture = (name) => new URL(`./fixtures/${name}`, import.meta.url).pathname;
const expected = JSON.parse(readFileSync(fixture('expected.json')));
await h5wasm.ready;

function open(name) {
  const file = new h5wasm.File(fixture(name), 'r');
  try {
    const info = describeFile(file);
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
  assert.equal(info.mask, '/MDHistoWorkspace/data/mask');
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
