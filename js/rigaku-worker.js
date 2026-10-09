// Module worker: reduces the frames of a Rigaku Oxford Diffraction (CrysAlisPro)
// experiment folder to an HKL volume (rigaku-reduce.js) and returns it as a Mantid
// MDHistoWorkspace file, which the page then opens like any other .nxs file.
//
// Messages in:  { type: 'scan', files: [{ file, path }] }
//               { type: 'scan-background', files } (air-scatter frames; files: [] clears them)
//               { type: 'reduce', options } (options.backgroundMode 'static' | 'rotation', backgroundBin
//               in degrees, backgroundScale: factor on the background)
// Messages out: { type: 'scanned', summary } | { type: 'background', summary } | { type: 'progress', label, fraction }
//               { type: 'done', blob, name, report } | { type: 'error', message }

import h5wasm from 'https://cdn.jsdelivr.net/npm/h5wasm@0.10.3/dist/esm/hdf5_hl.js';
import * as fmt from './rigaku-format.js';
import { cellFromUB } from './rigaku-geometry.js';
import { WorkerPool } from './rigaku-pool.js';
import { finishVolume, gridEdges, matchBackgroundRun, reduceRigaku, writeMantidMD } from './rigaku-reduce.js';

let exp = null, meta = null, bkg = null;

self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'scan') await scan(data.files);
    else if (data.type === 'scan-background') await scanBackground(data.files);
    else if (data.type === 'reduce') await reduce(data.options);
  } catch (err) {
    self.postMessage({ type: 'error', message: err?.message ?? String(err) });
  }
};

const text = async (f) => (f ? new TextDecoder('latin1').decode(await f.arrayBuffer()) : null);

async function scan(files) {
  bkg = null;
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

/**
 * Workers for the per-frame work (one per core, leaving one for this coordinator), or null:
 * then the reduction runs here, serially, with the same result.
 */
async function makePool() {
  const n = Math.min(12, (self.navigator?.hardwareConcurrency ?? 1) - 1);
  if (!(n >= 2) || typeof Worker === 'undefined') return null;
  let pool = null;
  try {
    pool = new WorkerPool(Array.from({ length: n }, () => new Worker(new URL('./rigaku-map-worker.js', import.meta.url), { type: 'module' })));
    await pool.ready();
    return pool;
  } catch {
    pool?.terminate();
    return null;
  }
}

/**
 * Air-scatter frames (crystal out of the beam, same setup): every frame of every run in the
 * folder. Their headers are compared with the experiment's here; reduceRigaku checks again.
 */
async function scanBackground(files) {
  if (!files.length) { bkg = null; self.postMessage({ type: 'background', summary: null }); return; }
  if (!exp) throw new Error('Choose the experiment folder first.');
  const b = fmt.scanExperiment(files.map(({ file }) => file));
  const list = [...b.runs].flatMap(([run, l]) => l.map((f) => ({ ...f, run })));
  if (!list.length) throw new Error('No frames (*_<run>_<frame>.rod_img) were found in the background folder.');
  const ref = fmt.parseRodHeader(await exp.runs.values().next().value[0].file.arrayBuffer());
  let exposure = 0;
  const differ = new Set(), runs = new Map();
  for (const f of list) {
    const h = fmt.parseRodHeader(await f.file.arrayBuffer());
    exposure += h.exposure;
    if (!runs.has(f.run)) runs.set(f.run, { run: f.run, axis: h.scanAxis, fixed: h.start.slice(0, 4), frames: [], start: h.start[h.scanAxis], end: h.end[h.scanAxis] });
    const r = runs.get(f.run);
    r.start = Math.min(r.start, h.start[r.axis]); r.end = Math.max(r.end, h.end[r.axis]);
    if (h.nx !== ref.nx || h.ny !== ref.ny) differ.add('detector size');
    if (Math.abs(h.distance - ref.distance) > 1e-3) differ.add(`distance (${h.distance} vs ${ref.distance} mm)`);
    if (Math.abs(h.start[1] - ref.start[1]) > 1e-3) differ.add(`2theta arm (${h.start[1]} vs ${ref.start[1]} deg)`);
    if (Math.abs(h.wavelengths.alpha12 - ref.wavelengths.alpha12) > 1e-6) differ.add('wavelength');
    if (h.pixelMM[0] !== ref.pixelMM[0] || h.binning[0] !== ref.binning[0] || h.binning[1] !== ref.binning[1]) differ.add('pixel size or binning');
  }
  bkg = differ.size ? null : { stem: b.stem, frames: list };
  // which experiment runs have a background run measured the same way (for a rotation-resolved background)
  const index = { runs: [...runs.values()] };
  const unmatched = [];
  for (const [run, l] of exp.runs) {
    const h = fmt.parseRodHeader(await l[0].file.arrayBuffer());
    if (!matchBackgroundRun(index, h.scanAxis, h.start)) unmatched.push(run);
  }
  self.postMessage({
    type: 'background',
    summary: {
      stem: b.stem, frames: list.length, runs: runs.size, exposure, matches: !differ.size, differences: [...differ], unmatchedRuns: unmatched,
      runList: index.runs.map((r) => ({ run: r.run, axis: r.axis, kappa: r.fixed[2], phi: r.fixed[3], omega: r.fixed[0], start: r.start, end: r.end })),
    },
  });
}

async function reduce(options) {
  if (!exp) throw new Error('Choose an experiment folder first.');
  const runs = new Map([...exp.runs].map(([r, list]) => [r, list.map(({ frame, file }) => ({ frame, read: () => file.arrayBuffer() }))]));
  const source = { runs, label: exp.stem, laue: meta.crystal.laue, monochromator: meta.par.monochromator, ubCandidates: meta.ubCandidates };
  const scale = options.backgroundScale ?? 1;
  if (bkg) {
    if (!(scale >= 0 && Number.isFinite(scale))) throw new Error('The background scale must be a number of 0 or more.');
    source.background = {
      frames: bkg.frames.map(({ run, frame, file }) => ({ run, frame, read: () => file.arrayBuffer() })),
      mode: options.backgroundMode ?? 'static', omegaBin: options.backgroundBin ?? 5,
    };
  }
  let last = 0;
  self.postMessage({ type: 'progress', label: 'Starting the workers', fraction: 0 });
  const pool = await makePool();
  let res;
  try {
    res = await reduceRigaku(source, fmt, options, (label, fraction) => {
      const now = performance.now();
      if (now - last > 150 || fraction >= 1) {
        last = now;
        self.postMessage({ type: 'progress', label, fraction: 0.97 * fraction });
      }
    }, pool);
  } finally {
    pool?.terminate();
  }
  self.postMessage({ type: 'progress', label: 'Writing the NeXus file', fraction: 0.97 });
  const { report, grid, model, ubMantid } = res;
  const v = finishVolume(res.acc, { inPlace: true, backgroundScale: scale });
  report.covered = v.covered;
  report.measuredZero = v.zeros;
  if (report.background) Object.assign(report.background, { stem: bkg.stem, scale, negativeVoxels: v.negative });
  const cell = report.grid.cell;
  const rows = [0, 3, 6].map((i) => report.grid.transform.slice(i, i + 3));
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
      xrd_cell_transform: `rows = output basis vectors in units of the refined cell (a b c): ${rows.map((r) => `(${r.join(', ')})`).join(' ')}; output indices = T x refined indices`,
      ...(report.background ? {
        xrd_background: `${report.background.mode === 'rotation' ? 'rotation-resolved' : 'static'} background from ${report.background.frames} frames `
          + `(${bkg.stem}, ${report.background.exposure} s), subtracted per pixel as counts before normalization, scale ${scale}: `
          + `${report.background.method}; ${report.background.errors}`,
      } : {}),
      xrd_reduction_report: JSON.stringify(report),
      xrd_source: `${exp.stem} (${report.frames} frames, runs ${report.runs.join(', ')})`,
    },
  });
  const blob = new Blob([FS.readFile(path)], { type: 'application/x-hdf5' });
  FS.unlink(path);
  const num = (x) => String(Number(x.toPrecision(6)));
  const cellTag = report.grid.multiplier ? `x${num(report.grid.multiplier)}` : `T${rows.map((r) => `(${r.map(num).join(',')})`).join('')}`;
  const bkgTag = report.background ? `_bkg${report.background.mode === 'rotation' ? 'rot' : ''}${scale === 1 ? '' : num(scale)}` : '';
  const name = `${exp.stem}_hkl_${cellTag}_step${Number(grid.step.toPrecision(6))}${bkgTag}.nxs`;
  self.postMessage({ type: 'done', blob, name, report });
}
