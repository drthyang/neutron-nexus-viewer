// The "Reindex with a new UB" dialog: load an ISAW UB file, show how its indexing relates to
// dataset A's, choose the new grid, resample A in its worker (worker.js, reindex.js), then
// open the resulting Mantid MDHistoWorkspace as a new dataset or to compare, or save it.

import { cellFromUB } from './nexus.js';
import { defaultSteps, defaultSubsamples, describeTransform, parseIsawUB, reindexGeometry, reindexGrid, reindexMatrix, reindexProblem } from './reindex.js';

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => (Number.isFinite(x) ? Number(x.toFixed(d)).toString().replace('-', '−') : '—');
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const cellText = (l) => `a ${fmt(l.a, 4)}, b ${fmt(l.b, 4)}, c ${fmt(l.c, 4)} Å, ${fmt(l.alpha, 2)}/${fmt(l.beta, 2)}/${fmt(l.gamma, 2)}°`;
const rows = (m, d) => [0, 3, 6].map((i) => `(${m.slice(i, i + 3).map((x) => fmt(x, d)).join(', ')})`).join(' ');
const AXES = ['H′', 'K′', 'L′'];
// 24 bytes per voxel while building (signal, σ², contributions), as for a 401³ volume.
const MAX_VOXELS = 6.5e7;

/**
 * `hooks`: datasetA() returns A's description (dims, lattice) or null; name() its file name;
 * post(message) sends to A's worker; openFile(file) opens a file as dataset A,
 * openCompare(file) as the next comparison dataset; canCompare() says whether one can be added.
 * Returns handle(message) for the worker's replies, refresh() for when A is ready and
 * reset() for when A is replaced.
 */
export function setupReindex(hooks) {
  const dlg = $('reindex');
  // ub: the loaded file, its UB and the transformation; box: the measured voxels in its indices.
  let ub = null, box = null, job = 0, busy = false, result = null;

  const note = (el, msg, kind = '') => {
    el.textContent = msg;
    el.className = `note ${kind}`;
    el.hidden = !msg;
  };
  const setBusy = (on) => {
    busy = on;
    $('ri-file').disabled = on;
    $('ri-cancel').disabled = $('ri-close').disabled = on;
    for (const el of $('ri-steps').querySelectorAll('input')) el.disabled = on;
    $('ri-run').disabled = on || !(currentGrid()?.voxels);
  };
  const progress = (label, fraction) => {
    $('ri-progress').hidden = false;
    $('ri-progress-label').textContent = label;
    $('ri-progress-pct').textContent = fraction > 0 ? `${Math.round(100 * fraction)}%` : '';
    $('ri-progress-fill').style.width = `${Math.round(100 * Math.min(1, Math.max(0, fraction)))}%`;
  };

  /** The panel card: what A is indexed with, or why it cannot be reindexed. */
  function refresh() {
    const A = hooks.datasetA(), why = reindexProblem(A);
    $('reindex-open').disabled = !!why;
    note($('reindex-status'), why || `Indexed with the file's UB: ${cellText(A.lattice)}.`);
  }

  // ---- UB file ----------------------------------------------------------------------
  $('reindex-open').onclick = $('ri-file').onclick = () => $('ri-input').click();
  $('ri-input').onchange = async () => {
    const file = $('ri-input').files[0];
    $('ri-input').value = '';
    if (!file) return;
    const A = hooks.datasetA(), why = reindexProblem(A);
    const where = dlg.open ? $('ri-note') : $('reindex-status');
    if (why) return note(where, why, 'error');
    try {
      load(file.name, parseIsawUB(await file.text()), A);
    } catch (err) {
      note(where, `${file.name}: ${err.message}`, 'error');
    }
  };

  function load(name, parsed, A) {
    const cellNew = cellFromUB(parsed.ub);
    if (!cellNew) throw new Error('its UB matrix gives no valid cell.');
    const T = reindexMatrix(A.lattice.ub, parsed.ub), geometry = reindexGeometry(A.dims, T);
    ub = { name, ...parsed, T, geometry, how: describeTransform(T, A.lattice.ub, parsed.ub), cellNew };
    box = null;
    result = null;
    const steps = defaultSteps(geometry.M, geometry.axes);
    stepInputs().forEach((el, i) => { el.value = String(steps[i]); });
    $('ri-nsub').value = defaultSubsamples(geometry.M, geometry.axes, steps);
    showSummary(A);
    $('ri-options').hidden = false;
    $('ri-result').hidden = true;
    $('ri-progress').hidden = true;
    note($('ri-note'), '');
    note($('ri-grid-note'), 'Finding the measured voxels in the new indices…');
    $('ri-run').disabled = true;
    hooks.post({ type: 'reindex-plan', id: ++job, T });
    if (!dlg.open) dlg.showModal();
  }

  function showSummary(A) {
    const { how, T, cellNew, cell: listed } = ub;
    const items = [
      ['Dataset A', `${esc(hooks.name())}: ${cellText(A.lattice)} (its UB)`],
      ['New UB', `${esc(ub.name)}: ${cellText(cellNew)}`],
    ];
    const off = listed && (['a', 'b', 'c'].some((k) => Math.abs(listed[k] / cellNew[k] - 1) > 0.005) || ['alpha', 'beta', 'gamma'].some((k) => Math.abs(listed[k] - cellNew[k]) > 0.5));
    if (off) items.push(['', `<span class="rk-bad">The file lists ${cellText(listed)}, which its UB does not give: is it an ISAW UB file?</span>`]);
    items.push(['Transformation', `<span class="ri-mono">H′ = T·H, T = ${rows(T, 4)}</span>`]);
    let relation;
    if (how.deviation < 0.1) {
      const identity = how.N.every((x, i) => x === (i % 4 ? 0 : 1));
      relation = `${identity ? 'T ≈ identity: the same setting' : `T ≈ <span class="ri-mono">${rows(how.N, 0)}</span>: another setting`}`
        + (how.angle !== null ? `, with the orientation turned by ${fmt(how.angle, 3)}°` : '');
    } else {
      relation = `<span class="rk-bad">T is not close to an integer matrix (off by up to ${fmt(how.deviation, 3)}): the new UB describes a differently oriented lattice. Both UBs must be in the same sample frame.</span>`;
    }
    if (how.det < 0) relation += ' <span class="rk-bad">det T &lt; 0: the new UB has the opposite handedness.</span>';
    items.push(['Relation', relation]);
    $('ri-summary').innerHTML = items.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
  }

  // ---- grid ---------------------------------------------------------------------------
  const stepInputs = () => [0, 1, 2].map((i) => $(`ri-step-${i}`));

  /** The grid of the current voxel sizes, an Error saying what is wrong, or null before the plan. */
  function currentGrid() {
    if (!ub || !box) return null;
    const steps = stepInputs().map((el) => Number(el.value)), nsub = Math.round(Number($('ri-nsub').value));
    if (!steps.every((s) => s > 0)) return new Error('Each voxel size must be a positive number.');
    if (!(nsub >= 1 && nsub <= 8)) return new Error('Sub-samples per axis: 1 to 8.');
    const grid = reindexGrid(box, steps);
    if (grid.voxels > MAX_VOXELS) return new Error(`${grid.shape.join(' × ')} = ${(grid.voxels / 1e6).toFixed(0)} M voxels is too many to build in the browser; use larger voxels.`);
    return Object.assign(grid, { nsub });
  }

  function showGrid() {
    const grid = currentGrid();
    if (!grid) return;
    if (grid instanceof Error) {
      note($('ri-grid-note'), grid.message, 'error');
      $('ri-run').disabled = true;
      return;
    }
    const ranges = AXES.map((x, i) => `${x} ${fmt(grid.edges[i][0], 3)} … ${fmt(grid.edges[i].at(-1), 3)}`).join(', ');
    const offCentre = AXES.filter((_, i) => Math.abs(1 / grid.steps[i] - Math.round(1 / grid.steps[i])) > 1e-6);
    const count = grid.voxels < 1e6 ? grid.voxels.toLocaleString() : `${(grid.voxels / 1e6).toFixed(2)} M`;
    note($('ri-grid-note'), `${grid.shape.join(' × ')} voxels (${count}), covering the measured voxels: ${ranges}. `
      + `${grid.nsub ** 3} sub-samples per voxel.${offCentre.length ? ` Integer ${offCentre.join(', ')} are not on voxel centres (1/voxel is not an integer).` : ''}`);
    $('ri-run').disabled = busy;
  }
  for (const el of $('ri-steps').querySelectorAll('input')) el.oninput = showGrid;

  // ---- run ----------------------------------------------------------------------------
  const stem = () => hooks.name().replace(/\.[^.]+$/, '');
  $('ri-run').onclick = () => {
    const grid = currentGrid();
    if (!grid || grid instanceof Error) return;
    const A = hooks.datasetA(), c = ub.cellNew;
    result = null;
    $('ri-result').hidden = true;
    note($('ri-note'), '');
    setBusy(true);
    progress('Starting', 0);
    hooks.post({
      type: 'reindex', id: ++job, T: ub.T, nsub: grid.nsub, ub: ub.ub, cell: [c.a, c.b, c.c, c.alpha, c.beta, c.gamma],
      grid: { min: grid.min, steps: grid.steps, shape: grid.shape, voxels: grid.voxels, edges: grid.edges },
      title: `${stem()}: reindexed with ${ub.name} in the NeXus Viewer`,
      logs: {
        reindex_source: hooks.name(),
        reindex_ub_file: ub.name,
        reindex_transform: `new indices = T x old indices, T rows: ${rows(ub.T, 8).replace(/−/g, '-')}`,
        reindex_previous_ub: `rows: ${rows(A.lattice.ub, 10).replace(/−/g, '-')}`,
        reindex_subsamples: grid.nsub,
      },
    });
  };

  function showResult(msg) {
    const file = new File([msg.blob], `${stem()}_reindexed.nxs`, { type: 'application/x-hdf5' });
    result = { file };
    const grid = currentGrid();
    $('ri-report').innerHTML = [
      ['File', `${esc(file.name)} (${(file.size / 1e6).toFixed(1)} MB)`],
      ['Grid', `${grid.shape.join(' × ')} voxels of ${grid.steps.map((s) => fmt(s, 4)).join(', ')} r.l.u.; ${msg.covered.toLocaleString()} measured (${(100 * msg.covered / msg.voxels).toFixed(1)} %)`],
      ['Errors', msg.errors ? 'propagated from the file’s uncertainties' : `none${msg.errorsNote ? ` (${esc(msg.errorsNote)})` : ' in the file'}`],
      ['Time', `${msg.seconds.toFixed(1)} s`],
    ].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('ri-open-compare').hidden = !hooks.canCompare();
    $('ri-result').hidden = false;
  }

  /** A reply from A's worker; replies to an earlier request are dropped. */
  function handle(msg) {
    if (msg.id !== job) return;
    if (msg.type === 'reindex-plan') {
      box = msg.box;
      showGrid();
    } else if (msg.type === 'progress-reindex') {
      progress(msg.label, msg.fraction);
    } else if (msg.type === 'reindex-file') {
      setBusy(false);
      progress('Done', 1);
      showResult(msg);
    } else if (msg.type === 'reindex-error') {
      setBusy(false);
      $('ri-progress').hidden = true;
      note($('ri-note'), msg.message, 'error');
    }
  }

  // ---- result and dialog plumbing --------------------------------------------------------
  $('ri-open-new').onclick = () => {
    if (!result) return;
    dlg.close();
    hooks.openFile(result.file);
  };
  $('ri-open-compare').onclick = () => {
    if (!result || !hooks.canCompare()) return;
    dlg.close();
    hooks.openCompare(result.file);
  };
  $('ri-download').onclick = () => {
    if (!result) return;
    const a = document.createElement('a');
    a.href = URL.createObjectURL(result.file);
    a.download = result.file.name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 10000);
  };
  $('ri-cancel').onclick = $('ri-close').onclick = () => { if (!busy) dlg.close(); };
  dlg.addEventListener('cancel', (e) => { if (busy) e.preventDefault(); });

  /** Dataset A was replaced: forget the UB and any reply still on its way. */
  function reset() {
    job++;
    ub = box = result = null;
    if (busy) setBusy(false);
    if (dlg.open) dlg.close();
    $('ri-options').hidden = $('ri-result').hidden = $('ri-progress').hidden = true;
  }

  return { handle, refresh, reset };
}
