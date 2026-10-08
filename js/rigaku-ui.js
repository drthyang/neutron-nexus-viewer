// The "Reduce Rigaku XRD" dialog: choose a CrysAlisPro experiment folder, pick runs and
// options, reduce it in a worker (rigaku-worker.js), then open the resulting Mantid
// MDHistoWorkspace as dataset A or as a comparison dataset, and offer it for download.

import { isHKL } from './nexus.js';
import { det3, parseFraction, parseTransform, transformCell } from './rigaku-geometry.js';

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => (Number.isFinite(x) ? Number(x.toFixed(d)).toString() : '—');
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// Output-cell presets: rows = new basis vectors in units of the refined cell's (a b c)
const PRESETS = {
  ub: [1, 0, 0, 0, 1, 0, 0, 0, 1],
  x2: [2, 0, 0, 0, 2, 0, 0, 0, 2],
  ortho: [1, 0, 0, 1, 2, 0, 0, 0, 1], // a' = a, b' = a + 2b, c' = c (hexagonal -> C-centred orthorhombic)
};
const isHexagonal = (c) => c && Math.abs(c[3] - 90) < 1 && Math.abs(c[4] - 90) < 1 && Math.abs(c[5] - 120) < 1 && Math.abs(c[0] - c[1]) / c[0] < 0.01;

/**
 * `hooks`: openFile(file) opens a file as dataset A; openCompare(file) as the next comparison
 * dataset; canCompare() says whether one can be added; datasetA() returns A's description
 * (dims with bin edges, lattice) or null.
 */
export function setupRigaku(hooks) {
  const dlg = $('rigaku');
  let worker = null, summary = null, result = null;

  const setBusy = (busy) => {
    $('rk-run').disabled = busy || !summary;
    $('rk-folder').disabled = busy;
    $('rk-cancel').textContent = busy ? 'Cancel' : 'Close';
  };
  const progress = (label, fraction) => {
    $('rk-progress').hidden = false;
    $('rk-progress-label').textContent = label;
    $('rk-progress-pct').textContent = fraction > 0 ? `${Math.round(100 * fraction)}%` : '';
    $('rk-progress-fill').style.width = `${Math.round(100 * Math.min(1, Math.max(0, fraction)))}%`;
  };
  const note = (msg, kind = '') => {
    $('rk-note').textContent = msg;
    $('rk-note').className = `note ${kind}`;
    $('rk-note').hidden = !msg;
  };
  const stop = () => { worker?.terminate(); worker = null; };

  const open = () => {
    const A = hooks.datasetA();
    const hkl = A && A.dims?.every(isHKL);
    $('rk-match-wrap').hidden = !hkl;
    $('rk-openas-wrap').hidden = !hooks.canCompare();
    if (!hkl && $('rk-cell').dataset.value === 'match') setPreset('ub');
    if (hooks.canCompare()) select($('rk-openas'), 'compare');
    if (summary) updateDefaults(); // dataset A may have changed (or finished loading) since
    dlg.showModal();
  };

  // ---- folder -------------------------------------------------------------------
  $('rk-folder').onclick = () => $('rk-input').click();
  $('rk-input').onchange = () => {
    const files = [...$('rk-input').files].map((file) => ({ file, path: file.webkitRelativePath || file.name }));
    $('rk-input').value = '';
    if (files.length) scan(files);
  };
  function scan(files) {
    stop();
    summary = null;
    result = null;
    $('rk-result').hidden = true;
    $('rk-options').hidden = true;
    note('');
    setBusy(true);
    progress(`Reading ${files.length.toLocaleString()} files`, 0);
    worker = new Worker(new URL('./rigaku-worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = ({ data }) => handlers[data.type]?.(data);
    worker.onerror = (e) => {
      e.preventDefault();
      fail(e.message || 'The reduction worker could not start (h5wasm is loaded from cdn.jsdelivr.net).');
    };
    worker.postMessage({ type: 'scan', files });
  }

  const fail = (message) => {
    setBusy(false);
    $('rk-progress').hidden = true;
    note(message, 'error');
  };

  const handlers = {
    error: ({ message }) => fail(message),
    progress: ({ label, fraction }) => progress(label, fraction),
    scanned: ({ summary: s }) => {
      summary = s;
      setBusy(false);
      $('rk-progress').hidden = true;
      showSummary(s);
    },
    done: ({ blob, name, report }) => {
      result = { file: new File([blob], name, { type: 'application/x-hdf5' }), report };
      setBusy(false);
      progress('Done', 1);
      showReport(report, name, blob.size);
      if ($('rk-autoopen').checked) openResult();
    },
  };

  function showSummary(s) {
    const T = s.temperature ? `${s.temperature[0]}–${s.temperature[1]} K` : '—';
    const cells = Object.entries(s.cells).map(([k, c]) => `${esc(k)}: a ${fmt(c[0])}, b ${fmt(c[1])}, c ${fmt(c[2])} Å, ${fmt(c[3], 2)}/${fmt(c[4], 2)}/${fmt(c[5], 2)}°`).join('<br>');
    const mono = s.monochromator?.theta ? `graphite, 2θ<sub>m</sub> ${fmt(2 * s.monochromator.theta, 2)}°, ${esc(s.monochromator.plane ?? '')}` : 'none recorded (unpolarized)';
    $('rk-summary').innerHTML = [
      ['Experiment', esc(s.stem)],
      ['Frames', `${s.frames.toLocaleString()} in ${s.runs.length} runs${s.preRuns ? ` (+${s.preRuns} pre-experiment runs, not used)` : ''}`],
      ['Detector', `${esc(s.detector)}, ${fmt(s.distance, 2)} mm, 2θ arm ${fmt(s.theta, 2)}°`],
      ['Wavelength', `${fmt(s.wavelength, 5)} Å`],
      ['Monochromator', mono],
      ['Temperature', T],
      ['Laue class', esc(s.laue ?? '—')],
      ['Cells (UBs)', cells || '—'],
      ['Decoding', s.decodeCheck ? 'first frame matches its header statistics' : '<span class="rk-bad">first frame does NOT match its header statistics</span>'],
    ].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('rk-runs').innerHTML = s.runs.map((r) => `<label class="rk-run"><input type="checkbox" data-run="${r.run}" checked>
      <span><b>${r.run}</b> ${r.frames} × ${fmt(r.width, 2)}° ${r.axis} ${fmt(r.start, 1)}→${fmt(r.end, 1)}°, κ ${fmt(r.kappa, 1)}°, φ ${fmt(r.phi, 2)}°, ${fmt(r.exposure, 2)} s${r.missing ? ` <span class="rk-bad">${r.missing} missing</span>` : ''}</span></label>`).join('');
    $('rk-options').hidden = false;
    const hex = isHexagonal(Object.values(s.cells)[0]);
    $('rk-ortho-wrap').hidden = !hex;
    if (!hex && $('rk-cell').dataset.value === 'ortho') setPreset('ub');
    updateDefaults();
  }

  // ---- output cell ----------------------------------------------------------------
  const matrixInputs = () => [...$('rk-matrix').querySelectorAll('input')].sort((a, b) => 3 * a.dataset.r + +a.dataset.c - (3 * b.dataset.r + +b.dataset.c));
  const setMatrix = (N) => matrixInputs().forEach((el, i) => { el.value = String(Number(N[i].toPrecision(6))); el.classList.remove('bad'); });
  const startCell = () => Object.values(summary?.cells ?? {})[0];
  function matchMultiple() {
    const A = hooks.datasetA(), c = startCell();
    return A?.lattice && c ? Math.max(1, Math.round(A.lattice.a / c[0])) : 1;
  }
  function setPreset(value) {
    select($('rk-cell'), value);
    if (value === 'match') { const n = matchMultiple(); setMatrix([n, 0, 0, 0, n, 0, 0, 0, n]); } else if (PRESETS[value]) setMatrix(PRESETS[value]);
  }
  /** The transformation in the matrix editor, or an Error describing what is wrong with it. */
  function readMatrix() {
    const els = matrixInputs();
    try {
      const N = parseTransform(els.map((el) => el.value));
      els.forEach((el) => el.classList.remove('bad'));
      return N;
    } catch (err) {
      els.forEach((el) => el.classList.toggle('bad', !Number.isFinite(parseFraction(el.value))));
      return err;
    }
  }

  function updateDefaults() {
    const A = hooks.datasetA();
    const mode = $('rk-cell').dataset.value;
    if (mode === 'match') setPreset('match'); // A may have changed since the preset was chosen
    const N = readMatrix(), c = startCell();
    if (N instanceof Error) {
      $('rk-cell-note').textContent = N.message;
      $('rk-cell-note').className = 'note error';
      return;
    }
    $('rk-cell-note').className = 'note';
    const d = det3(N), out = c ? transformCell(c, N) : null;
    const cellText = out ? `a′ ${fmt(out[0])}, b′ ${fmt(out[1])}, c′ ${fmt(out[2])} Å, ${fmt(out[3], 2)}/${fmt(out[4], 2)}/${fmt(out[5], 2)}°, ${fmt(d, 4)} × the cell volume` : '';
    if (mode === 'match' && A?.lattice) {
      const e = A.dims[0].edges;
      $('rk-step').value = fmt(e[1] - e[0], 4);
      $('rk-cell-note').textContent = `Output cell ${cellText} (dataset A: a ≈ ${fmt(A.lattice.a)} Å); voxel centres on A's.`;
    } else {
      if (!$('rk-step').dataset.touched) {
        const r = Math.cbrt(d), f = r >= 1 ? Math.round(r) : 1 / Math.max(1, Math.round(1 / r));
        $('rk-step').value = fmt(0.05 * f, 4);
      }
      $('rk-cell-note').textContent = out ? `Output cell ${cellText} (from the CrysAlis UB; refined in the reduction). Indices: H′ = T·(h, k, l).` : '';
    }
  }

  function options() {
    const runs = [...$('rk-runs').querySelectorAll('input[data-run]')].filter((i) => i.checked).map((i) => Number(i.dataset.run));
    if (!runs.length) throw new Error('Select at least one run.');
    const opts = {
      runs, step: Number($('rk-step').value), nSub: Math.max(1, Math.round(Number($('rk-nsub').value) || 5)),
      normalization: $('rk-norm').dataset.value, refine: $('rk-refine').checked,
    };
    if (!(opts.step > 0)) throw new Error('The voxel size must be positive.');
    const mode = $('rk-cell').dataset.value, A = hooks.datasetA();
    if (mode === 'match') setPreset('match');
    const N = readMatrix();
    if (N instanceof Error) throw new Error(`Output cell: ${N.message}`);
    opts.transform = N;
    if (mode === 'match') {
      if (!A?.lattice || !A.dims?.every(isHKL)) throw new Error('Dataset A is not open yet, or has no HKL axes and cell to match.');
      // A's edges may carry float32 rounding (0.10000038): keep 6 significant digits
      const round6 = (x) => Number(x.toPrecision(6));
      opts.origin = A.dims.map((d) => round6(0.5 * (d.edges[0] + d.edges[1])));
      if (!$('rk-step').dataset.touched) opts.step = round6(A.dims[0].edges[1] - A.dims[0].edges[0]);
    }
    return opts;
  }

  $('rk-run').onclick = () => {
    let opts;
    try { opts = options(); } catch (err) { note(err.message, 'error'); return; }
    note('');
    $('rk-result').hidden = true;
    setBusy(true);
    progress('Starting', 0);
    worker.postMessage({ type: 'reduce', options: opts });
  };

  function showReport(r, name, size) {
    const L = r.geometry?.levels?.L10;
    const rows = [
      ['File', `${esc(name)} (${(size / 1e6).toFixed(1)} MB)`],
      ['Frames', `${r.frames.toLocaleString()} (runs ${r.runs.join(', ')}); all match their header statistics; max ${r.frameCheck.maxPixelCounts.toLocaleString()} counts/pixel`],
      ['Mask', `${(100 * r.mask.maskedFraction).toFixed(1)} % of pixels (beamstop ${r.mask.beamstopPixels.toLocaleString()}, chip lines ${r.mask.boundaryColumns.length} + ${r.mask.boundaryRows.length}, dead ${r.mask.deadPixels})`],
      ['Peaks', `${r.peaks.found} found, ${r.peaks.clean} clear of the mask; ${r.indexed} indexed with ${esc(r.startUB)} (${Object.entries(r.ubCandidates).map(([k, v]) => `${esc(k)} ${(100 * v).toFixed(0)} %`).join(', ')})`],
      ['Geometry', `${L ? `refined: rms ${fmt(L.rmsX, 2)} / ${fmt(L.rmsY, 2)} px, ${fmt(L.rmsAngle, 3)}° (median ${fmt(L.medX, 2)} / ${fmt(L.medY, 2)} px, ${fmt(L.medAngle, 3)}°)${r.geometry.levels.L10a.driftRuns?.length ? `; drift in run ${r.geometry.levels.L10a.driftRuns.join(', ')}` : ''}` : 'header model, not refined'}; detector distance ${fmt(r.model.geometry.distance, 3)} mm (header calibration)`],
      ['Cell', `a ${fmt(r.cell[0], 4)}, b ${fmt(r.cell[1], 4)}, c ${fmt(r.cell[2], 4)} Å, ${fmt(r.cell[3], 2)}/${fmt(r.cell[4], 2)}/${fmt(r.cell[5], 2)}° (${esc(r.geometry.system ?? '')})`],
      ['Output cell', `a′ ${fmt(r.grid.cell[0], 4)}, b′ ${fmt(r.grid.cell[1], 4)}, c′ ${fmt(r.grid.cell[2], 4)} Å, ${fmt(r.grid.cell[3], 2)}/${fmt(r.grid.cell[4], 2)}/${fmt(r.grid.cell[5], 2)}°; ${r.grid.multiplier ? `${r.grid.multiplier} × the cell` : `T = ${[0, 3, 6].map((i) => `(${r.grid.transform.slice(i, i + 3).map((v) => fmt(v, 4)).join(', ')})`).join(' ')}`}`],
      ['Grid', `${r.grid.shape.join(' × ')} voxels of ${r.grid.step} r.l.u.; ${r.covered.toLocaleString()} measured (${r.measuredZero.toLocaleString()} with zero counts)`],
      ['Signal', esc(r.normalization)],
      ['Time', `${r.seconds.toFixed(0)} s`],
    ];
    $('rk-report').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('rk-result').hidden = false;
  }

  const download = (blob, name) => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 10000);
  };
  $('rk-download').onclick = () => result && download(result.file, result.file.name);
  $('rk-download-report').onclick = () => result && download(new Blob([JSON.stringify(result.report, null, 1)], { type: 'application/json' }),
    result.file.name.replace(/\.nxs$/, '_report.json'));
  $('rk-open').onclick = openResult;

  function openResult() {
    if (!result) return;
    const asCompare = !$('rk-openas-wrap').hidden && $('rk-openas').dataset.value === 'compare' && hooks.canCompare();
    dlg.close();
    if (asCompare) hooks.openCompare(result.file);
    else hooks.openFile(result.file);
  }

  // ---- dialog plumbing ------------------------------------------------------------------
  function select(seg, value) {
    seg.dataset.value = value;
    for (const b of seg.querySelectorAll('button')) b.classList.toggle('on', b.dataset.value === value);
  }
  for (const seg of dlg.querySelectorAll('.segmented')) {
    for (const b of seg.querySelectorAll('button')) {
      b.onclick = () => {
        if (seg.id === 'rk-cell') { setPreset(b.dataset.value); updateDefaults(); } else select(seg, b.dataset.value);
      };
    }
  }
  // editing the matrix leaves the presets (and their step defaults stay unless touched)
  for (const el of matrixInputs()) el.oninput = () => { select($('rk-cell'), 'custom'); updateDefaults(); };
  $('rk-step').oninput = () => { $('rk-step').dataset.touched = '1'; };
  $('rk-cancel').onclick = () => {
    if (worker && $('rk-cancel').textContent === 'Cancel') {
      stop();
      summary = null;
      $('rk-options').hidden = true;
      $('rk-progress').hidden = true;
      setBusy(false);
      note('Cancelled. Choose the folder again to restart.');
      return;
    }
    dlg.close();
  };
  $('rk-close').onclick = () => $('rk-cancel').click();
  dlg.addEventListener('cancel', (e) => { if (worker && $('rk-cancel').textContent === 'Cancel') e.preventDefault(); });
  for (const id of ['rigaku-open', 'rigaku-intro']) $(id)?.addEventListener('click', open);
  return { open };
}
