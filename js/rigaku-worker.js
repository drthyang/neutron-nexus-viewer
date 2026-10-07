// Module worker: reduces the frames of a Rigaku Oxford Diffraction (CrysAlisPro)
// experiment folder to an HKL volume (rigaku-reduce.js) and returns it as a Mantid
// MDHistoWorkspace file, which the page then opens like any other .nxs file.
//
// Messages in:  { type: 'scan', files: [{ file, path }] }
//               { type: 'reduce', options }
// Messages out: { type: 'scanned', summary } | { type: 'progress', label, fraction }
//               { type: 'done', blob, name, report } | { type: 'error', message }

import h5wasm from 'https://cdn.jsdelivr.net/npm/h5wasm@0.10.3/dist/esm/hdf5_hl.js';
import * as fmt from './rigaku-format.js';
import { cellFromUB } from './rigaku-geometry.js';
import { finishVolume, gridEdges, reduceRigaku, writeMantidMD } from './rigaku-reduce.js';

let exp = null, meta = null;

self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'scan') await scan(data.files);
    else if (data.type === 'reduce') await reduce(data.options);
  } catch (err) {
    self.postMessage({ type: 'error', message: err?.message ?? String(err) });
  }
};

const text = async (f) => (f ? new TextDecoder('latin1').decode(await f.arrayBuffer()) : null);

async function scan(files) {
  exp = fmt.scanExperiment(files.map(({ file }) => file));
  if (!exp.runs.size) throw new Error('No CrysAlisPro frames (*_<run>_<frame>.rod_img) were found in this folder.');
  const [crystal, par, cracker, coll] = await Promise.all([
    text(exp.meta.crystalIni).then((t) => (t ? fmt.parseCrystalIni(t) : {})),
    text(exp.meta.par).then((t) => (t ? fmt.parsePar(t) : {})),
    text(exp.meta.crackerPar).then((t) => (t ? fmt.parsePar(t) : {})),
    text(exp.meta.datacollIni).then((t) => (t ? fmt.parseDatacoll(t) : {})),
  ]);
  meta = {
    crystal, par, cracker, coll,
    ubCandidates: { 'crystal.ini Gral UB': crystal.gralUB, 'crystal.ini Mother of all UBs': crystal.motherUB, '_cracker.par': cracker.ub, '.par': par.ub },
  };
  const runs = [];
  let first = null;
  for (const [run, list] of exp.runs) {
    const head = fmt.parseRodHeader(await list[0].file.arrayBuffer());
    const last = fmt.parseRodHeader(await list.at(-1).file.arrayBuffer());
    first ??= head;
    const ax = head.scanAxis;
    runs.push({
      run, frames: list.length, missing: list.at(-1).frame - list[0].frame + 1 - list.length,
      axis: ax === 0 ? 'omega' : ax === 3 ? 'phi' : `axis ${ax}`, start: head.start[ax], end: last.end[ax],
      width: head.end[ax] - head.start[ax], exposure: head.exposure,
      omega: head.start[0], theta: head.start[1], kappa: head.start[2], phi: head.start[3],
    });
  }
  const buf = await exp.runs.values().next().value[0].file.arrayBuffer();
  const img = fmt.decodeTY6(buf, first);
  const lambda = first.wavelengths.alpha12;
  const cells = Object.fromEntries(Object.entries(meta.ubCandidates).filter(([, ub]) => ub).map(([k, ub]) => [k, cellFromUB(ub, lambda)]));
  self.postMessage({
    type: 'scanned',
    summary: {
      stem: exp.stem, runs, preRuns: exp.pre.size, frames: runs.reduce((s, r) => s + r.frames, 0),
      detector: `${first.nx} × ${first.ny} px, ${first.pixelMM[0]} mm`, distance: first.distance, theta: first.start[1],
      wavelength: lambda, laue: meta.crystal.laue, temperature: meta.coll.temperature, collected: meta.coll.start,
      monochromator: meta.par.monochromator, cells, decodeCheck: fmt.checkStats(first, img), metaFiles: Object.keys(exp.meta),
    },
  });
}

async function reduce(options) {
  if (!exp) throw new Error('Choose an experiment folder first.');
  const runs = new Map([...exp.runs].map(([r, list]) => [r, list.map(({ frame, file }) => ({ frame, read: () => file.arrayBuffer() }))]));
  const source = { runs, label: exp.stem, laue: meta.crystal.laue, monochromator: meta.par.monochromator, ubCandidates: meta.ubCandidates };
  let last = 0;
  const res = await reduceRigaku(source, fmt, options, (label, fraction) => {
    const now = performance.now();
    if (now - last > 150 || fraction >= 1) {
      last = now;
      self.postMessage({ type: 'progress', label, fraction: 0.97 * fraction });
    }
  });
  self.postMessage({ type: 'progress', label: 'Writing the NeXus file', fraction: 0.97 });
  const { report, grid, model, ubMantid } = res;
  const v = finishVolume(res.acc, { inPlace: true });
  report.covered = v.covered;
  report.measuredZero = v.zeros;
  const m = report.grid.multiplier;
  const cell = [m * model.cell[0], m * model.cell[1], m * model.cell[2], model.cell[3], model.cell[4], model.cell[5]];
  const { FS } = await h5wasm.ready;
  const path = '/rigaku-reduced.nxs';
  const T = meta.coll.temperature;
  writeMantidMD(h5wasm, path, {
    shape: [grid.shape[2], grid.shape[1], grid.shape[0]], edges: gridEdges(grid), signal: v.signal, errors2: v.errors2, events: v.events,
    ub: ubMantid, cell, title: `${exp.stem}: Rigaku XRD reduced in the NeXus Viewer`,
    logs: {
      wavelength: { value: model.g.wavelength, units: 'Angstrom' },
      ...(T ? { temperature_min: { value: T[0], units: 'K' }, temperature_max: { value: T[1], units: 'K' } } : {}),
      xrd_signal_definition: report.normalization,
      xrd_reduction_report: JSON.stringify(report),
      xrd_source: `${exp.stem} (${report.frames} frames, runs ${report.runs.join(', ')})`,
    },
  });
  const blob = new Blob([FS.readFile(path)], { type: 'application/x-hdf5' });
  FS.unlink(path);
  const name = `${exp.stem}_hkl_x${m}_step${grid.step}.nxs`;
  self.postMessage({ type: 'done', blob, name, report });
}
