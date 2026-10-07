import assert from 'node:assert/strict';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { describeFile, loadVariance, loadVolume } from '../js/nexus.js';
import * as fmt from '../js/rigaku-format.js';
import {
  M_CRYSALIS, bMatrix, cellFromUB, crystalSystem, gonio, inv3, levenbergMarquardt, mul, mulv, pixelX, predict, prepare, rot, transpose,
} from '../js/rigaku-geometry.js';
import {
  PeakSearch, accumulateFrame, detectorMask, finishVolume, gridEdges, makeAccumulators, makeGrid, reduceRigaku, writeMantidMD,
} from '../js/rigaku-reduce.js';
import { encodeTY6, rng, simulateExperiment, writeFrame } from './rigaku-synth.js';

await h5wasm.ready;

test('TY6 round trip: every encoding branch decodes to the original pixels', () => {
  const nx = 53, ny = 7, img = new Int32Array(nx * ny), r = rng(3);
  for (let i = 0; i < img.length; i++) {
    const u = r();
    img[i] = u < 0.5 ? Math.floor(r() * 3) : u < 0.8 ? Math.floor(r() * 200) : u < 0.95 ? Math.floor(r() * 30000) : Math.floor(r() * 3e6);
  }
  img[0] = 5e6; img[nx] = -40; img[2 * nx] = 200; // first-pixel long, plain and short forms
  for (let x = 0; x < nx; x++) img[3 * nx + x] = 7; // a row of zero differences (nbit 0)
  const header = { nx, ny, start: [10, 20, 54, 30], end: [10.5, 20, 54, 30], zeroCorr: [90, 0.05, 0, 0], exposure: 3, pixelMM: 0.1,
    origin: [26.5, 3.25], distance: 45.28, detRot: [-0.44, 0.02, 0], alpha: 90, wavelengths: [0.7093, 0.71359, 0.71073, 0.63229] };
  const bytes = writeFrame(img, header);
  const h = fmt.parseRodHeader(bytes);
  assert.equal(h.compression, 'TY6');
  assert.deepEqual([h.nx, h.ny, h.nheader], [nx, ny, 6576]);
  assert.equal(h.scanAxis, 0);
  assert.ok(Math.abs(h.start[0] - 10) < 1e-9 && Math.abs(h.end[0] - 10.5) < 1e-9 && Math.abs(h.zeroCorr[0] - 90) < 1e-9);
  assert.deepEqual([h.distance, h.exposure, h.alpha, h.wavelengths.alpha12], [45.28, 3, 90, 0.71073]);
  assert.deepEqual(h.origin, [26.5, 3.25]);
  const out = fmt.decodeTY6(bytes, h);
  assert.deepEqual(Array.from(out), Array.from(img));
  assert.ok(fmt.checkStats(h, out));
  // a corrupted pixel no longer matches the header statistics
  out[5] += 1;
  assert.ok(!fmt.checkStats(h, out));
});

test('experiment folder: frames by run, pre-experiment runs apart, metadata by name', () => {
  const names = ['X_1_2.rod_img', 'X_1_1.rod_img', 'X_2_1.rod_img', 'pre_X_1_1.rod_img', 'X.par', 'X_cracker.par', 'pre_X.par', 'X_crystal.ini',
    'pre_X_crystal.ini', 'X_datacoll.ini', 'notes.txt'];
  const exp = fmt.scanExperiment(names.map((name) => ({ name })));
  assert.equal(exp.stem, 'X');
  assert.deepEqual([...exp.runs.keys()], [1, 2]);
  assert.deepEqual(exp.runs.get(1).map((f) => f.frame), [1, 2]);
  assert.equal(exp.pre.size, 1);
  assert.deepEqual(Object.fromEntries(Object.entries(exp.meta).map(([k, f]) => [k, f.name])),
    { par: 'X.par', crackerPar: 'X_cracker.par', crystalIni: 'X_crystal.ini', datacollIni: 'X_datacoll.ini' });
  assert.throws(() => fmt.scanExperiment([{ name: 'A_1_1.rod_img' }, { name: 'B_1_1.rod_img' }]), /several experiments/);
});

test('CrysAlis text files: UB, monochromator, Laue class, instrument model', () => {
  const par = '³XCALIBUR SYSTEM\n³   - ALPHA (DEG)   90.00000 BETA (DEG)    0.00000\n'
    + '³   - MONOCHROMATOR DVALUE (ANG)    3.35400 MONOCHROMATOR THETA (DEG)    6.06977\n'
    + 'CRYSTALLOGRAPHY UB   1.0E-002   0.0   0.0   0.0   2.0E-002   0.0   0.0   0.0   3.0E-002\n'
    + 'ROTATION MONOCHROMATOR   3.354000 E1E3PLANE  0.50\nEXPERIMENT TEMPERATURE    294.000\n';
  const p = fmt.parsePar(par);
  assert.deepEqual(p.ub, [0.01, 0, 0, 0, 0.02, 0, 0, 0, 0.03]);
  assert.deepEqual(p.monochromator, { dvalue: 3.354, theta: 6.06977, plane: 'E1E3PLANE' });
  assert.deepEqual([p.temperature, p.alpha, p.beta], [294, 90, 0]);
  const doubles = new Float64Array(32);
  doubles[24] = -0.44; doubles[25] = 0.02; doubles[27] = 389.8; doubles[28] = 192.6; doubles[31] = 45.49;
  const hex = Array.from(new Uint8Array(doubles.buffer), (b) => b.toString(16).padStart(2, '0')).join('');
  const ini = `[Gral UB]\r\nmatrix=0.1 0 0 0 0.1 0 0 0 0.05 0\r\n[Symmetry]\r\nlaue class="6/m"\r\n[Mother of all lattices]\r\nInstrumentModel binary_0=${hex}\r\n`;
  const c = fmt.parseCrystalIni(ini);
  assert.deepEqual(c.gralUB, [0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.05]);
  assert.equal(c.laue, '6/m');
  assert.deepEqual(c.motherModel, { origin: [389.8, 192.6], distance: 45.49, detRot: [-0.44, 0.02] });
  assert.deepEqual(fmt.parseDatacoll('[Sample T in K]\r\nSample T in K min max=294.090000 295.120000\r\n').temperature, [294.09, 295.12]);
});

test('geometry: cells, prediction and its inverse agree', () => {
  const cell = [4.0, 4.0, 5.0, 90, 90, 120];
  assert.ok(cellFromUB(bMatrix(...cell)).every((v, i) => Math.abs(v - cell[i]) < 1e-9));
  assert.equal(crystalSystem('6/m', cell), 'hexagonal');
  assert.equal(crystalSystem('6/m', [4, 4.2, 5, 90, 90, 120]), 'triclinic'); // a != b: not that setting
  const g = { wavelength: 0.71073, ox: 388.8, oy: 192.4, distance: 45.28, pixelMM: 0.1, d1: -0.44, d2: 0.02, theta: 20, thetaOffset: 0.05,
    alpha: 90, omegaOffset: 90, kappaOffset: 0.02, phiOffset: 0, tx: -0.05, ty: 0.06, tz: 0.01 };
  const G = prepare(g);
  const U = mul(rot([1, 0, 0], 12), rot([0, 0, 1], 33));
  const ub = mul(U, bMatrix(...cell)).map((v) => v * g.wavelength);
  const inv = mul(inv3(ub), M_CRYSALIS);
  let n = 0;
  for (const hkl of [[1, 0, 0], [2, -1, 1], [1, 1, 2], [0, 2, -3], [3, -2, 0]]) {
    for (const near of [-60, 0, 60]) {
      const p = predict(G, mulv(mul(transpose(M_CRYSALIS), ub), hkl), [0, 20, 54, 30], 0, near);
      if (!p) continue;
      const R = gonio(G, p.angle, 54, 30);
      const h = mulv(inv, mulv(transpose(R), pixelX(G, p.i, p.j, R)));
      h.forEach((v, i) => assert.ok(Math.abs(v - hkl[i]) < 1e-9, `${hkl} at ${p.angle}: ${h}`));
      n++;
    }
  }
  assert.ok(n >= 5);
});

test('Levenberg-Marquardt recovers a model and down-weights outliers', () => {
  const xs = Array.from({ length: 40 }, (_, i) => i / 4);
  const ys = xs.map((x) => 2.5 * Math.exp(-0.3 * x) + 0.7);
  ys[7] += 50; // one wild point
  const fit = levenbergMarquardt((p) => Float64Array.from(xs, (x, i) => p[0] * Math.exp(-p[1] * x) + p[2] - ys[i]), [1, 0.1, 0],
    { steps: [1e-6, 1e-6, 1e-6], fScale: 0.1 });
  assert.ok(Math.abs(fit.p[0] - 2.5) < 0.02 && Math.abs(fit.p[1] - 0.3) < 0.01 && Math.abs(fit.p[2] - 0.7) < 0.02, String(fit.p));
});

test('gridding: shares landing in one voxel are summed before squaring', () => {
  const grid = makeGrid([-1, -1, -1], [1, 1, 1], 0.5);
  assert.deepEqual(grid.shape, [5, 5, 5]);
  const acc = makeAccumulators(grid);
  const I = [1, 0, 0, 0, 1, 0, 0, 0, 1];
  // pixel 0: all five sub-samples in the voxel at the origin; pixel 1: split 2/3 between two voxels
  const x = Float64Array.of(0, 0, 0, 0.2, 0, 0);
  const A = [I, I, I.map((v, i) => (i === 0 ? 1.5 : v)), I.map((v, i) => (i === 0 ? 1.5 : v)), I.map((v, i) => (i === 0 ? 1.5 : v))];
  accumulateFrame(acc, grid, x, Float64Array.of(10, 10), Float64Array.of(3, 3), A);
  const at = (h, k, l) => ((l + 2) * 5 + (k + 2)) * 5 + (h + 2);
  assert.equal(acc.S[at(0, 0, 0)], 10 + 10 * 2 / 5);
  assert.ok(Math.abs(acc.E2[at(0, 0, 0)] - (10 + 10 * 4 / 25)) < 1e-5);
  assert.equal(acc.S[at(1, 0, 0)], 10 * 3 / 5); // 0.2 * 1.5 = 0.3 -> voxel H = 0.5
  assert.ok(Math.abs(acc.E2[at(1, 0, 0)] - 10 * 9 / 25) < 1e-5);
  assert.equal(acc.W[at(0, 0, 0)], 3 + 3 * 2 / 5);
  assert.equal(acc.N[at(0, 0, 0)], 2);
  const v = finishVolume(acc);
  assert.ok(Number.isNaN(v.signal[at(-1, -1, -1)]) && v.events[at(-1, -1, -1)] === 0);
  assert.equal(v.covered, 2);
  assert.ok(Math.abs(v.signal[at(0, 0, 0)] - 14 / 4.2) < 1e-12);
});

test('detector mask: beamstop with penumbra, chip-boundary triplet, dead pixel, border', () => {
  const nx = 200, ny = 120, sum = new Float64Array(nx * ny).fill(1000), counted = new Uint8Array(nx * ny).fill(1);
  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      const r = Math.hypot(x - 60, y - 60);
      if (r < 6) sum[y * nx + x] = 10;
      else if (r < 10) sum[y * nx + x] = 500;
    }
  }
  for (let y = 0; y < ny; y++) { sum[y * nx + 140] *= 0.6; sum[y * nx + 139] *= 1.3; sum[y * nx + 141] *= 1.3; }
  counted[30 * nx + 30] = 0;
  const { mask, info } = detectorMask(sum, counted, nx, ny);
  assert.ok(mask[60 * nx + 60] && mask[60 * nx + 68], 'umbra and penumbra');
  assert.ok(!mask[60 * nx + 75], 'clear pixels outside the penumbra');
  assert.deepEqual(info.boundaryColumns, [139, 140, 141]);
  assert.ok(mask[30 * nx + 30] && info.deadPixels === 1);
  assert.ok(mask[0] && mask[(ny - 1) * nx + 5] && mask[10 * nx + nx - 1]);
});

test('peak search: one spot over three frames, net-count weighted centroid', () => {
  const nx = 40, ny = 30, mask = new Uint8Array(nx * ny), bg = new Float64Array(nx * ny).fill(1);
  const ps = new PeakSearch(nx, ny, mask, bg);
  const weights = [0.25, 1, 0.5];
  for (const w of weights) {
    const img = new Int32Array(nx * ny).fill(1);
    for (const [dx, dy, a] of [[0, 0, 400], [1, 0, 200], [0, 1, 200], [-1, 0, 200], [0, -1, 200]]) img[(12 + dy) * nx + 20 + dx] = 1 + Math.round(a * w);
    ps.addFrame(img);
  }
  const peaks = ps.peaks();
  assert.equal(peaks.length, 1);
  const p = peaks[0];
  assert.ok(Math.abs(p.x - 20) < 1e-9 && Math.abs(p.y - 12) < 1e-9);
  const zc = weights.reduce((s, w, k) => s + w * (k + 0.5), 0) / weights.reduce((s, w) => s + w, 0);
  assert.ok(Math.abs(p.z - zc) < 1e-3, `${p.z} vs ${zc}`);
  assert.equal(p.frames, 3);
});

test('Mantid MDHistoWorkspace writer: the viewer reads dims, NaN coverage, uncertainties and the cell', () => {
  const grid = makeGrid([-1, -2, -3], [1, 2, 3], 1); // 3 x 5 x 7
  const n = 3 * 5 * 7, signal = new Float64Array(n), errors2 = new Float64Array(n), events = new Float64Array(n);
  for (let v = 0; v < n; v++) { signal[v] = v % 4 === 0 ? NaN : v / 10; errors2[v] = v % 4 === 0 ? NaN : 0.01 * v; events[v] = v % 4 === 0 ? 0 : 3; }
  const ub = mul(rot([0, 0, 1], 10), bMatrix(8, 8, 10, 90, 90, 120));
  const path = '/tmp/rigaku-writer-test.nxs';
  writeMantidMD(h5wasm, path, { shape: [7, 5, 3], edges: gridEdges(grid), signal, errors2, events, ub, cell: [8, 8, 10, 90, 90, 120],
    logs: { wavelength: { value: 0.71073, units: 'Angstrom' }, note: 'test' }, title: 'writer test' });
  const f = new h5wasm.File(path, 'r');
  try {
    const info = describeFile(f);
    assert.deepEqual(info.shape, [7, 5, 3]);
    assert.deepEqual(info.dims.map((d) => d.label), ['H', 'K', 'L']);
    assert.ok(info.dims.every((d) => d.frame === 'HKL'));
    assert.deepEqual(Array.from(info.dims[2].edges), [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]);
    assert.ok(Math.abs(info.lattice.a - 8) < 1e-9 && Math.abs(info.lattice.c - 10) < 1e-9 && Math.abs(info.lattice.gamma - 120) < 1e-9);
    const { volume, stats } = loadVolume(f, info);
    assert.equal(stats.valid, n - Math.ceil(n / 4));
    assert.ok(Number.isNaN(volume[0]) && Math.abs(volume[1] - 0.1) < 1e-6);
    const variance = loadVariance(f, info);
    assert.ok(Math.abs(variance[5] - 0.05) < 1e-7);
    assert.equal(f.get('MDHistoWorkspace').attrs.QConvention.value, 'Crystallography');
    assert.equal(f.get('MDHistoWorkspace/data/signal').attrs.axes.value, 'D2:D1:D0');
    assert.equal(f.get('MDHistoWorkspace/experiment0/logs/note/value').value[0], 'test');
  } finally {
    f.close();
  }
});

test('end to end on a simulated experiment: geometry refined, Bragg peaks at integer HKL', async () => {
  const sim = simulateExperiment();
  // CrysAlis's matrix is slightly off: 0.4 deg misoriented and a 0.5 % larger cell
  const lambda = sim.model.g.wavelength;
  const U = mul(sim.model.ub.map((v) => v / lambda), inv3(bMatrix(...sim.model.cell)));
  const ubStart = mul(mul(rot([0.3, -0.5, 0.8], 0.4), U), bMatrix(4.02, 4.02, 5.025, 90, 90, 120)).map((v) => v * lambda);
  const runs = new Map();
  for (const fr of sim.frames) {
    if (!runs.has(fr.run)) runs.set(fr.run, []);
    runs.get(fr.run).push({ frame: fr.frame, read: async () => fr.bytes.buffer.slice(0) });
  }
  const res = await reduceRigaku({ runs, ubCandidates: { start: ubStart }, laue: '6/m', label: 'simulated' }, fmt,
    { multiplier: 1, step: 0.05, nSub: 3, normalization: 'rate' });
  const r = res.report;
  assert.ok(r.peaks.clean >= 40, `peaks: ${JSON.stringify(r.peaks)}`);
  assert.ok(r.indexed / r.peaks.clean > 0.9);
  const L = r.geometry.levels.L10;
  assert.ok(L.rmsX < 0.3 && L.rmsY < 0.3 && L.rmsAngle < 0.1, JSON.stringify(L));
  assert.ok(Math.abs(res.model.cell[0] / 4.0 - 1) < 0.005 && Math.abs(res.model.cell[2] / 5.0 - 1) < 0.005, String(res.model.cell));
  // the strongest voxels sit on integer HKL
  const v = finishVolume(res.acc), [nH, nK] = res.grid.shape, step = res.grid.step, min = res.grid.min;
  const order = Array.from(v.signal.keys()).filter((i) => v.signal[i] === v.signal[i]).sort((a, b) => v.signal[b] - v.signal[a]).slice(0, 15);
  for (const i of order) {
    const ih = i % nH, ik = Math.floor(i / nH) % nK, il = Math.floor(i / (nH * nK));
    const hkl = [min[0] + ih * step, min[1] + ik * step, min[2] + il * step];
    assert.ok(hkl.every((c) => Math.abs(c - Math.round(c)) <= step + 1e-9), `bright voxel at ${hkl}`);
  }
});
