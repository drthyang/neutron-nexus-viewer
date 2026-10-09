import assert from 'node:assert/strict';
import { unlinkSync } from 'node:fs';
import { tmpdir } from 'node:os';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { cellFromUB, describeFile, loadVariance, loadVolume } from '../js/nexus.js';
import {
  defaultSteps, defaultSubsamples, describeTransform, measuredBox, parseIsawUB, reindexGeometry, reindexGrid, reindexMatrix,
  reindexProblem, reindexVolume,
} from '../js/reindex.js';
import { inv3, mul, rot } from '../js/rigaku-geometry.js';
import { writeMantidMD } from '../js/rigaku-reduce.js';

const fixture = (name) => new URL(`./fixtures/${name}`, import.meta.url).pathname;
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

/** An ISAW UB file as Mantid's SaveIsawUB writes it: row j holds column j of the ISAW-frame UB (row r = Mantid row (r + 2) mod 3). */
function isawText(ub, cell = null) {
  const rows = [0, 1, 2].map((j) => [0, 1, 2].map((r) => ub[3 * ((r + 2) % 3) + j].toFixed(8)).join(' '));
  if (cell) rows.push(`${cell.join(' ')} 100.0`, '0 0 0 0 0 0 0');
  return `${rows.join('\n')}\n\n\nThe above matrix is the Transpose of the UB Matrix.\n`;
}

const close = (a, b, tol, label) => a.forEach((x, i) => assert.ok(Math.abs(x - b[i]) <= tol, `${label}[${i}]: ${x} vs ${b[i]}`));

/** Every new voxel against the old voxel holding its centre h′ (old index h = Tinv·h′ by integer arithmetic on the bins). */
function checkNearest(src, grid, out, oldOf, label) {
  const { info, volume, variance } = src;
  const [nH, nK, nL] = grid.shape;
  let compared = 0;
  for (let c = 0, o = 0; c < nL; c++) {
    for (let b = 0; b < nK; b++) {
      for (let a = 0; a < nH; a++, o++) {
        const h = [grid.min[0] + a * grid.steps[0], grid.min[1] + b * grid.steps[1], grid.min[2] + c * grid.steps[2]];
        const u = oldOf(h); // old display coordinates of the centre
        const j = info.dims.map((d, k) => Math.round((u[k] - d.edges[0]) / (d.edges[1] - d.edges[0]) - 0.5));
        const inside = j.every((x, k) => x >= 0 && x < info.dims[k].edges.length - 1);
        const v = inside ? volume[j[0] + info.shape[2] * (j[1] + info.shape[1] * j[2])] : NaN;
        if (Number.isNaN(v)) {
          assert.ok(Number.isNaN(out.signal[o]), `${label} voxel ${o} should be empty`);
        } else {
          assert.equal(out.signal[o], v, `${label} voxel ${o}`);
          if (variance) assert.ok(Math.abs(out.errors2[o] - variance[j[0] + info.shape[2] * (j[1] + info.shape[1] * j[2])]) < 1e-12, `${label} σ² ${o}`);
          compared++;
        }
      }
    }
  }
  return compared;
}

test('ISAW UB files: the frame swap and transpose of Mantid LoadIsawUB', () => {
  // Orthorhombic, a* along Mantid x, b* along y (up), c* along z (the beam). In ISAW's frame
  // (x the beam, z up) c* is along x, b* along z and a* along y, and the file is UBᵀ.
  const [a, b, c] = [4, 5, 8];
  const text = `0 ${1 / a} 0\n0 0 ${1 / b}\n${1 / c} 0 0\n${a} ${b} ${c} 90 90 90 160\n0 0 0 0 0 0 0\n\nThe above matrix is the Transpose of the UB Matrix.\n`;
  const { ub, cell } = parseIsawUB(text);
  close(ub, [1 / a, 0, 0, 0, 1 / b, 0, 0, 0, 1 / c], 1e-15, 'diagonal UB');
  assert.deepEqual(cell, { a, b, c, alpha: 90, beta: 90, gamma: 90 });

  // A general UB survives the SaveIsawUB layout, and its cell is the one the file lists.
  const { info } = open('mdhisto_small.nxs');
  const l = info.lattice;
  const back = parseIsawUB(isawText(l.ub, [l.a, l.b, l.c, l.alpha, l.beta, l.gamma]));
  close(back.ub, l.ub, 1e-8, 'round trip');
  const fromUB = cellFromUB(back.ub);
  for (const k of ['a', 'b', 'c', 'alpha', 'beta', 'gamma']) assert.ok(Math.abs(fromUB[k] - back.cell[k]) < 1e-5, k);

  assert.throws(() => parseIsawUB('hello\n1 2 3'), /not an ISAW UB file/);
  assert.throws(() => parseIsawUB('1 0 0\n0 1 0\n0 0 0\n'), /singular/);
});

test('the transformation between two UBs, and the misorientation once in the same setting', () => {
  const ub = open('mdhisto_small.nxs').info.lattice.ub;
  const N = [1, 0, 0, 1, 2, 0, 0, 0, 1]; // orthohexagonal: h′ = h, k′ = h + 2k, l′ = l
  const T = reindexMatrix(ub, mul(ub, inv3(N)));
  close(T, N, 1e-12, 'T');
  const same = describeTransform(T, ub, mul(ub, inv3(N)));
  assert.deepEqual(same.N, N);
  assert.ok(same.deviation < 1e-12 && same.angle < 1e-6 && Math.abs(same.det - 2) < 1e-12);

  const turned = mul(rot([0, 1, 0], 0.3), mul(ub, inv3(N)));
  const d = describeTransform(reindexMatrix(ub, turned), ub, turned);
  assert.deepEqual(d.N, N);
  assert.ok(d.deviation > 1e-4 && d.deviation < 0.05, `deviation ${d.deviation}`);
  assert.ok(Math.abs(d.angle - 0.3) < 1e-6, `angle ${d.angle}`);
});

test('the same UB reproduces the volume voxel by voxel, with σ² of a voxel sampled several times unchanged', () => {
  const src = open('mdhisto_small.nxs'), { info, volume, variance } = src;
  assert.equal(reindexProblem(info), '');
  const T = reindexMatrix(info.lattice.ub, info.lattice.ub);
  const { TW, M, axes } = reindexGeometry(info.dims, T);
  const steps = defaultSteps(M, axes);
  close(steps, [0.1, 0.1, 0.25], 1e-12, 'steps');
  assert.equal(defaultSubsamples(M, axes, steps), 2);
  const grid = reindexGrid(measuredBox(volume, info.shape, axes, TW), steps);
  // The measured voxels' box lies on the old bins.
  grid.edges.forEach((e, i) => e.forEach((x) => assert.ok(info.dims[i].edges.some((y) => Math.abs(x - y) < 1e-9), `edge ${x}`)));
  for (const nsub of [1, 2, 3]) {
    const out = reindexVolume(volume, variance, info.shape, axes, M, grid, nsub);
    const n = checkNearest(src, grid, out, (h) => h, `nsub ${nsub}`);
    assert.equal(n, info.shape.reduce((p, x) => p * x, 1) * src.stats.fraction);
    assert.equal(out.covered, n);
  }
  // The same UB read back from an ISAW file (8 decimals) gives the same grid.
  const read = reindexGeometry(info.dims, reindexMatrix(info.lattice.ub, parseIsawUB(isawText(info.lattice.ub)).ub));
  const again = reindexGrid(measuredBox(volume, info.shape, read.axes, read.TW), defaultSteps(read.M, read.axes));
  assert.deepEqual(again.shape, grid.shape);
  again.edges.forEach((e, i) => close(e, grid.edges[i], 1e-12, `edges ${i}`));
});

test('a 90° relabelling (h′ = k, k′ = −h) moves every voxel to its new indices', () => {
  const src = open('mdhisto_small.nxs'), { info, volume, variance } = src;
  const N = [0, 1, 0, -1, 0, 0, 0, 0, 1];
  const ubNew = mul(info.lattice.ub, inv3(N));
  const T = reindexMatrix(info.lattice.ub, ubNew);
  close(T, N, 1e-12, 'T');
  const { TW, M, axes } = reindexGeometry(info.dims, T);
  const steps = defaultSteps(M, axes);
  close(steps, [0.1, 0.1, 0.25], 1e-12, 'steps');
  const grid = reindexGrid(measuredBox(volume, info.shape, axes, TW), steps);
  const out = reindexVolume(volume, variance, info.shape, axes, M, grid, 2);
  const n = checkNearest(src, grid, out, ([h, k, l]) => [-k, h, l], 'relabelled');
  assert.equal(n, src.stats.valid);
});

test('axes stored in another order (a NEBULA3D volume, L first) come out as H, K, L', () => {
  const src = open('nebula3d_small.h5'), { info, volume, variance } = src;
  assert.deepEqual(info.dims.map((d) => d.label), ['L', 'K', 'H']);
  const { TW, M, axes } = reindexGeometry(info.dims, reindexMatrix(info.lattice.ub, info.lattice.ub));
  const steps = defaultSteps(M, axes);
  close(steps, [0.2, 0.1, 0.25], 1e-12, 'steps');
  const grid = reindexGrid(measuredBox(volume, info.shape, axes, TW), steps);
  const out = reindexVolume(volume, variance, info.shape, axes, M, grid, 1);
  assert.equal(checkNearest(src, grid, out, ([h, k, l]) => [l, k, h], 'nebula3d'), src.stats.valid);
});

test('a coarser grid averages the old voxels it holds, with σ² = Σσ²/N² over the measured ones', () => {
  const n = 8, step = 0.25, e = Array.from({ length: n + 1 }, (_, k) => -1 + k * step);
  const dims = ['H', 'K', 'L'].map((x, i) => ({ label: x, edges: e, frame: 'HKL', units: 'r.l.u.', basis: { vec: [0, 1, 2].map((j) => +(i === j)), letter: x } }));
  const shape = [n, n, n], volume = new Float32Array(n ** 3), variance = new Float32Array(n ** 3);
  let seed = 7;
  const rand = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  for (let v = 0; v < volume.length; v++) {
    volume[v] = rand() < 0.15 ? NaN : Math.round(1000 * rand()) / 8;
    variance[v] = Math.round(100 * rand()) / 16;
  }
  const { TW, M, axes } = reindexGeometry(dims, [1, 0, 0, 0, 1, 0, 0, 0, 1]);
  const grid = reindexGrid(measuredBox(volume, shape, axes, TW), [0.5, 0.5, 0.5]);
  assert.deepEqual(grid.shape, [5, 5, 5]); // centres −1 … 1: the outer ones hold one old layer, the others two
  const out = reindexVolume(volume, variance, shape, axes, M, grid, 2);
  for (let c = 0, o = 0; c < 5; c++) {
    for (let b = 0; b < 5; b++) {
      for (let a = 0; a < 5; a++, o++) {
        let sum = 0, e2 = 0, N = 0;
        const span = (k) => [2 * k - 1, 2 * k].filter((j) => j >= 0 && j < n);
        for (const z of span(c)) for (const y of span(b)) for (const x of span(a)) {
          const v = volume[x + n * (y + n * z)];
          if (Number.isNaN(v)) continue;
          sum += v;
          e2 += variance[x + n * (y + n * z)];
          N++;
        }
        if (!N) { assert.ok(Number.isNaN(out.signal[o]), `voxel ${o} empty`); continue; }
        assert.ok(Math.abs(out.signal[o] - sum / N) < 1e-12, `voxel ${o}: ${out.signal[o]} vs ${sum / N}`);
        assert.ok(Math.abs(out.errors2[o] - e2 / (N * N)) < 1e-12, `σ² ${o}`);
        assert.equal(out.events[o], N);
      }
    }
  }
});

test('a reindexed volume written in the Mantid layout reads back with its grid, UB and values', () => {
  const src = open('mdhisto_small.nxs'), { info, volume, variance } = src;
  const N = [1, 0, 0, 1, 2, 0, 0, 0, 1], ubNew = mul(rot([0, 0, 1], 0.2), mul(info.lattice.ub, inv3(N)));
  const T = reindexMatrix(info.lattice.ub, ubNew);
  const { TW, M, axes } = reindexGeometry(info.dims, T);
  const steps = defaultSteps(M, axes);
  close(steps, [0.1, 0.2, 0.25], 1e-12, 'steps'); // k′ = h + 2k: one old K bin is 0.2 in k′
  const grid = reindexGrid(measuredBox(volume, info.shape, axes, TW), steps);
  const out = reindexVolume(volume, variance, info.shape, axes, M, grid, defaultSubsamples(M, axes, steps));
  assert.ok(out.covered > 0);
  const c = cellFromUB(ubNew), path = `${tmpdir()}/reindex-test-${process.pid}.nxs`;
  writeMantidMD(h5wasm, path, {
    shape: [grid.shape[2], grid.shape[1], grid.shape[0]], edges: grid.edges, signal: out.signal, errors2: out.errors2, events: out.events,
    ub: ubNew, cell: [c.a, c.b, c.c, c.alpha, c.beta, c.gamma], logs: { reindex_transform: T.join(' ') },
  });
  try {
    const back = (() => {
      const file = new h5wasm.File(path, 'r');
      try {
        const i = describeFile(file);
        return { info: i, ...loadVolume(file, i), variance: loadVariance(file, i) };
      } finally {
        file.close();
      }
    })();
    assert.deepEqual(back.info.dims.map((d) => d.label), ['H', 'K', 'L']);
    back.info.dims.forEach((d, i) => close(d.edges, grid.edges[i], 1e-12, `edges ${i}`));
    close(back.info.lattice.ub, ubNew, 1e-15, 'UB');
    assert.equal(back.stats.valid, out.covered);
    out.signal.forEach((v, k) => {
      if (Number.isNaN(v)) assert.ok(Number.isNaN(back.volume[k]), `voxel ${k}`);
      else {
        assert.equal(back.volume[k], Math.fround(v), `voxel ${k}`);
        assert.equal(back.variance[k], Math.fround(out.errors2[k]), `σ² ${k}`);
      }
    });
  } finally {
    unlinkSync(path);
  }
});

test('datasets that cannot be reindexed say why', () => {
  assert.match(reindexProblem(null), /Open a dataset/);
  assert.match(reindexProblem(open('nxdata_small.h5').info), /r\.l\.u\./);
  const { info } = open('mdhisto_small.nxs');
  assert.match(reindexProblem({ ...info, lattice: { ...info.lattice, ub: undefined } }), /no UB matrix/);
  const edges = [...info.dims[0].edges];
  edges[3] += 0.01;
  assert.match(reindexProblem({ ...info, dims: [{ ...info.dims[0], edges }, ...info.dims.slice(1)] }), /uniform/);
});
