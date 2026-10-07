// The "Reduce Rigaku XRD" dialog: choose a CrysAlisPro experiment folder, pick runs and
// options, reduce it in a worker (rigaku-worker.js), then open the resulting Mantid
// MDHistoWorkspace as dataset A or as a comparison dataset, and offer it for download.

import { isHKL } from './nexus.js';

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => (Number.isFinite(x) ? Number(x.toFixed(d)).toString() : '—');
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

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
    if (!hkl && $('rk-cell').dataset.value === 'match') select($('rk-cell'), 'ub');
    if (hooks.canCompare()) select($('rk-openas'), 'compare');
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
    updateDefaults();
  }

  function updateDefaults() {
    const A = hooks.datasetA();
    const mode = $('rk-cell').dataset.value;
    const firstCell = Object.values(summary?.cells ?? {})[0];
    if (mode === 'match' && A?.lattice) {
      const e = A.dims[0].edges;
      $('rk-step').value = fmt(e[1] - e[0], 4);
      $('rk-cell-note').textContent = `Output cell a ≈ ${fmt(A.lattice.a)} Å (dataset A); voxel centres on A's.`;
    } else {
      const m = mode === 'x2' ? 2 : 1;
      if (!$('rk-step').dataset.touched) $('rk-step').value = fmt(0.05 * m, 3);
      $('rk-cell-note').textContent = firstCell ? `Output cell a ≈ ${fmt(m * firstCell[0])} Å, c ≈ ${fmt(m * firstCell[2])} Å (CrysAlis UB; refined in the reduction).` : '';
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
    if (mode === 'x2') opts.multiplier = 2;
    else if (mode === 'ub') opts.multiplier = 1;
    else if (mode === 'match' && A?.lattice) {
      opts.targetA = A.lattice.a;
      opts.origin = A.dims.map((d) => 0.5 * (d.edges[0] + d.edges[1]));
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
      ['Geometry', L ? `refined: rms ${fmt(L.rmsX, 2)} / ${fmt(L.rmsY, 2)} px, ${fmt(L.rmsAngle, 3)}° (median ${fmt(L.medX, 2)} / ${fmt(L.medY, 2)} px, ${fmt(L.medAngle, 3)}°)${r.geometry.levels.L10a.driftRuns?.length ? `; drift in run ${r.geometry.levels.L10a.driftRuns.join(', ')}` : ''}` : 'header model, not refined'],
      ['Cell', `a ${fmt(r.cell[0], 4)}, b ${fmt(r.cell[1], 4)}, c ${fmt(r.cell[2], 4)} Å, ${fmt(r.cell[3], 2)}/${fmt(r.cell[4], 2)}/${fmt(r.cell[5], 2)}° (${esc(r.geometry.system ?? '')}); output = ${r.grid.multiplier} × this cell`],
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
    for (const b of seg.querySelectorAll('button')) b.onclick = () => { select(seg, b.dataset.value); if (seg.id === 'rk-cell') updateDefaults(); };
  }
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
