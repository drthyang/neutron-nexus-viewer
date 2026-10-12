import assert from 'node:assert/strict';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

import h5wasm from 'h5wasm/node';

import { callTool, TOOLS, toolDefinitions, Workspace } from '../agent/index.js';
import { describeFile, loadVolume } from '../js/nexus.js';

const fixture = (name) => new URL(`./fixtures/${name}`, import.meta.url).pathname;
const demo = new URL('../examples/demo_10K.nxs', import.meta.url).pathname;
const expected = JSON.parse(readFileSync(fixture('expected.json')));
const tmp = mkdtempSync(join(tmpdir(), 'nexus-agent-'));

async function call(ws, name, args) {
  const r = await callTool(name, args, { workspace: ws });
  assert.ok(r.ok, `${name}: ${r.error}`);
  return r.value;
}

/** The float32 values of a .npy file the tools wrote. */
function readNpy(path) {
  const bytes = readFileSync(path), len = bytes.readUInt16LE(8);
  return new Float32Array(bytes.buffer.slice(bytes.byteOffset + 10 + len, bytes.byteOffset + bytes.length));
}

test('open_volume describes a Mantid file and gives the grid NEXPLAN takes', async () => {
  const ws = new Workspace();
  const v = await call(ws, 'open_volume', { path: fixture('mdhisto_small.nxs') });
  assert.equal(v.volume_id, 'mdhisto_small');
  assert.equal(v.format, 'Mantid MDHistoWorkspace (SaveMD)');
  assert.deepEqual(v.shape, [7, 9, 11]);
  assert.deepEqual(v.axes.map((a) => a.label), ['H', 'K', 'L']);
  assert.deepEqual(v.planes.map((p) => p.plane), ['HK', 'HL', 'KL']);
  assert.equal(v.cell.source, 'UB');
  assert.deepEqual(v.nexplan_grid.axes, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
  assert.deepEqual(v.nexplan_grid.binning, ['-0.55,0.1,0.55', '-0.45,0.1,0.45', '-0.875,0.25,0.875']);
  assert.match(v.uncertainties, /^variances/);
  // A second open of the same file gets its own id; describe_volume lists both.
  await call(ws, 'open_volume', { path: fixture('mdhisto_small.nxs') });
  const list = await call(ws, 'describe_volume', {});
  assert.deepEqual(list.volumes.map((x) => x.volume_id), ['mdhisto_small', 'mdhisto_small_2']);
});

test('open_volume reads NEBULA3D and NXdata files, and refuses others', async () => {
  const ws = new Workspace();
  assert.equal((await call(ws, 'open_volume', { path: fixture('nebula3d_small.h5') })).format, 'NEBULA3D volume');
  const dpdf = await call(ws, 'open_volume', { path: fixture('nebula3d_dpdf_small.h5') });
  assert.equal(dpdf.format, 'NEBULA3D 3D-ΔPDF');
  assert.match(dpdf.nexplan_grid.unavailable, /hkl/);
  assert.equal((await call(ws, 'open_volume', { path: fixture('nxdata_small.h5') })).format, 'NXdata');
  const bad = await callTool('open_volume', { path: fixture('expected.json') }, { workspace: ws });
  assert.equal(bad.ok, false);
  assert.match(bad.error, /not an HDF5/);
});

test('slice_volume gives the slices the Python reference gives', async () => {
  const ws = new Workspace();
  await call(ws, 'open_volume', { path: fixture('mdhisto_small.nxs'), id: 'v' });
  for (const c of expected.mdhisto_small.slice(0, 6)) {
    await call(ws, 'set_symmetry', { volume_id: 'v', ...(c.inversion ? { laue_class: '-1' } : { operations: 'none' }) });
    const out = join(tmp, `slice-${c.fixed}-${c.center}-${c.inversion}.npy`);
    const s = await call(ws, 'slice_volume', { volume_id: 'v', fixed_axis: c.fixed, center: c.center, thickness: c.thickness, image: false, output_path: out });
    assert.deepEqual(s.shape, c.shape);
    assert.equal(s.slab.bins, c.ids.length);
    const values = readNpy(out);
    c.values.forEach((v, i) => {
      if (v === null) assert.ok(Number.isNaN(values[i]));
      else assert.ok(Math.abs(values[i] - v) <= 1e-5 * Math.max(1, Math.abs(v)), `pixel ${i}: ${values[i]} vs ${v}`);
    });
  }
});

test('slice_volume returns a PNG in true geometry, and the viewer\'s Auto color range', async () => {
  const ws = new Workspace();
  await call(ws, 'open_volume', { path: demo, id: 'cold' });
  await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '6/mmm' });
  const r = await callTool('slice_volume', { volume_id: 'cold', fixed_axis: 'L' }, { workspace: ws });
  assert.ok(r.ok);
  assert.equal(r.images.length, 1);
  assert.deepEqual([...r.images[0].data.slice(1, 4)].map((b) => String.fromCharCode(b)).join(''), 'PNG');
  assert.equal(r.value.plane, 'HK');
  assert.ok(Math.abs(r.value.angle_between_axes_deg - 60) < 1e-6);
  assert.equal(r.value.color.min, 0);
  assert.equal(r.value.color.colormap, 'magma');
  assert.ok(r.value.image.width > r.value.image.height, 'a 60° parallelogram is wider than tall');
});

test('symmetry, mask and their errors behave as in the viewer', async () => {
  const ws = new Workspace(), ctx = { workspace: ws };
  await call(ws, 'open_volume', { path: demo, id: 'cold' });
  let r = await callTool('build_mask', { volume_id: 'cold', outlier_k: 5 }, ctx);
  assert.match(r.error, /at least 3 operations/);
  r = await callTool('slice_volume', { volume_id: 'cold', fixed_axis: 'L', removed_only: true }, ctx);
  assert.match(r.error, /needs a mask/);
  const sym = await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '6/mmm' });
  assert.equal(sym.order, 24);
  assert.ok(sym.metric_change < 1e-9);
  const wrong = await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '4/mmm' });
  assert.match(wrong.warning, /not symmetries of this lattice/);
  r = await callTool('set_symmetry', { volume_id: 'cold', laue_class: '432' }, ctx);
  assert.ok(r.details.presets.includes('m-3m'));
  await call(ws, 'set_symmetry', { volume_id: 'cold', operations: 'k,-h-k,l; -h,-k,-l; -k,-h,-l' });
  assert.equal((await call(ws, 'describe_volume', { volume_id: 'cold' })).symmetry.group, '-31m');
  const mask = await call(ws, 'build_mask', { volume_id: 'cold', edge_radius: 1, outlier_k: 5, output_path: join(tmp, 'mask.npy') });
  assert.ok(mask.coverage_edge_voxels > 0 && mask.outlier_voxels > 0);
  const bytes = readFileSync(mask.written.path);
  assert.equal(bytes.length - 10 - bytes.readUInt16LE(8), 101 ** 3);
  r = await callTool('build_mask', { volume_id: 'cold', edge_radius: 1, output_path: join(tmp, 'mask.npy') }, ctx);
  assert.match(r.error, /exists/);
  r = await callTool('slice_volume', { volume_id: 'cold', fixed_axis: 'L', bogus: 1 }, ctx);
  assert.match(r.error, /Invalid arguments/);
  r = await callTool('slice_volume', { volume_id: 'warm', fixed_axis: 'L' }, ctx);
  assert.deepEqual(r.details.open, ['cold']);
});

test('read_points gives the value a one-bin slice shows at that pixel', async () => {
  const ws = new Workspace();
  await call(ws, 'open_volume', { path: demo, id: 'cold' });
  await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '6/mmm' });
  const out = join(tmp, 'hk0.npy');
  await call(ws, 'slice_volume', { volume_id: 'cold', fixed_axis: 'L', center: 0, image: false, output_path: out });
  const slice = readNpy(out);
  const p = await call(ws, 'read_points', { volume_id: 'cold', points: [[1, 0, 0], [0.5, 0.5, 0], [9, 0, 0]], box: 0 });
  // Pixel (row = K bin, column = H bin) of H = 1, K = 0 and H = K = 0.5 on the -5.05..5.05 grid in steps of 0.1.
  assert.ok(Math.abs(p.points[0].value - slice[50 * 101 + 60]) < 1e-3);
  assert.ok(Math.abs(p.points[1].value - slice[55 * 101 + 55]) < 1e-3);
  assert.equal(p.points[2].voxel, null);
  assert.equal(p.summary.points.with_data, 2);
});

test('read_points takes a NEXPLAN bragg_positions file and checks its grid', async () => {
  const ws = new Workspace(), ctx = { workspace: ws };
  await call(ws, 'open_volume', { path: demo, id: 'cold' });
  await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '6/mmm' });
  // P6₃/mmc forbids hhl with odd l.
  const nodes = { schema: 'nexplan-bragg-nodes/1', grid: { axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]] }, nodes: [
    { hkl: [1, 0, 0], grid: [1, 0, 0], class: 'present' },
    { hkl: [1, 1, 1], grid: [1, 1, 1], class: 'systematic absence' },
  ] };
  writeFileSync(join(tmp, 'nodes.json'), JSON.stringify(nodes));
  const r = await call(ws, 'read_points', { volume_id: 'cold', points_path: join(tmp, 'nodes.json') });
  assert.deepEqual(Object.keys(r.summary), ['present', 'systematic absence']);
  assert.ok(r.summary.present.median_value > 10 * r.summary['systematic absence'].median_value);
  writeFileSync(join(tmp, 'other.json'), JSON.stringify({ ...nodes, grid: { axes: [[1, 1, 0], [-1, 1, 0], [0, 0, 1]] } }));
  const bad = await callTool('read_points', { volume_id: 'cold', points_path: join(tmp, 'other.json') }, ctx);
  assert.match(bad.error, /computed for axes/);
});

test('line_cut and powder_average give profiles with the processing applied', async () => {
  const ws = new Workspace();
  await call(ws, 'open_volume', { path: fixture('mdhisto_small.nxs'), id: 'v' });
  const cut = await call(ws, 'line_cut', { volume_id: 'v', start: [0, 0, 0.5], end: [0, 0, -0.5], output_path: join(tmp, 'cut.txt') });
  assert.match(cut.note, /swapped/);
  assert.equal(cut.axis, 'L');
  assert.equal(cut.width.unit, 'Å⁻¹');
  assert.ok(cut.points.some((p) => Number.isFinite(p.sigma)), 'σ from errors_squared');
  assert.match(readFileSync(cut.written.path, 'utf8'), /^# Written by NeXus Viewer agent tools/);
  const iq = await call(ws, 'powder_average', { volume_id: 'v', bins: '0.05', output_path: join(tmp, 'iq.csv') });
  assert.equal(iq.q_unit, 'Å⁻¹, with 2π');
  assert.ok(iq.shells.every((s) => s.coverage >= 0 && Number.isFinite(s.coverage)));
  assert.ok(iq.shells.length > 3);
  const bad = await callTool('powder_average', { volume_id: 'v', bins: '0' }, { workspace: ws });
  assert.match(bad.error, /cannot be 0/);
});

test('export_nebula3d writes the file NEBULA3D reads, with the processing recorded', async () => {
  const ws = new Workspace();
  await call(ws, 'open_volume', { path: demo, id: 'cold' });
  await call(ws, 'set_symmetry', { volume_id: 'cold', laue_class: '6/mmm' });
  const out = join(tmp, 'cold_sym.nxs');
  const e = await call(ws, 'export_nebula3d', { volume_id: 'cold', output_path: out });
  assert.deepEqual(e.shape_hkl, [101, 101, 101]);
  assert.equal(e.attributes.symmetry, '6/mmm');
  await h5wasm.ready;
  const f = new h5wasm.File(out, 'r');
  try {
    const info = describeFile(f);
    assert.equal(info.signal, '/entry/data');
    assert.equal(f.get('entry').attrs.symmetry.value, '6/mmm');
    assert.ok(loadVolume(f, info).stats.valid > 0);
  } finally {
    f.close();
  }
  const again = await callTool('export_nebula3d', { volume_id: 'cold', output_path: out }, { workspace: ws });
  assert.match(again.error, /exists/);
  // The example has no uncertainties: no sigma, and NEBULA3D's √|data| is said.
  assert.equal(e.sigma, false);
  assert.match(e.attributes.uncertainties, /none written/);
  // A Mantid file's errors_squared become /entry/sigma, σ of each orbit's mean.
  await call(ws, 'open_volume', { path: fixture('mdhisto_small.nxs'), id: 'v' });
  const withSigma = await call(ws, 'export_nebula3d', { volume_id: 'v', output_path: join(tmp, 'v.nxs') });
  assert.equal(withSigma.sigma, true);
  assert.match(withSigma.attributes.uncertainties, /^from errors_squared/);
  const g = new h5wasm.File(withSigma.written.path, 'r');
  try {
    const entry = g.get('entry');
    assert.deepEqual(entry.get('sigma').shape, entry.get('data').shape);
    const sigma = entry.get('sigma').value, mask = entry.get('mask').value;
    assert.ok(sigma.some((x, i) => mask[i] && x > 0) && sigma.every((x, i) => mask[i] || x === 0));
  } finally {
    g.close();
  }
});

test('every tool has a JSON Schema and is documented in docs/AGENT.md', () => {
  const docs = readFileSync(new URL('../docs/AGENT.md', import.meta.url), 'utf8');
  const defs = toolDefinitions();
  assert.equal(defs.length, TOOLS.length);
  for (const d of defs) {
    assert.equal(d.input_schema.type, 'object', d.name);
    assert.equal(d.input_schema.additionalProperties, false, `${d.name} refuses unknown arguments`);
    assert.ok(docs.includes(`\`${d.name}\``), `${d.name} is missing from docs/AGENT.md`);
  }
});
