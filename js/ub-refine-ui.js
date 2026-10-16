// The "Refine the UB" dialog: find the Bragg peaks of dataset A in its worker (worker.js,
// ub-refine.js), fit an orientation and a cell constrained to a crystal system here, and
// hand the refined UB to the reindexing dialog, or save it as an ISAW UB file.

import { formatIsawUB, reindexProblem } from './reindex.js';
import { UB_SYSTEMS, formatSU, indexShift, peaksCSV, refineUB, suggestSystem, voxelQ } from './ub-refine.js';
import { cellFromUB, inv3, mul } from './rigaku-geometry.js';

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => (Number.isFinite(x) ? Number(x.toFixed(d)).toString().replace('-', '−') : '—');
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const KEYS = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'];

/** A cell with its standard uncertainties: lengths in Å, then angles; values held fixed by the system without decimals. */
function cellText(cell, su = null) {
  const v = (i, d) => (su && su[i] > 0 ? formatSU(cell[i], su[i], d) : su ? String(Number(cell[i].toFixed(d))) : cell[i].toFixed(d));
  return `a ${v(0, 4)}, b ${v(1, 4)}, c ${v(2, 4)} Å; α ${v(3, 3)}, β ${v(4, 3)}, γ ${v(5, 3)}°`.replace(/-/g, '−');
}

/**
 * `hooks`: datasetA() returns A's description (dims, lattice) or null; name() its file name;
 * post(message) sends to A's worker; laue() the Laue class applied (a preset's name) or
 * null; reindexWith(name, parsed) opens the reindexing dialog with a UB ({ ub, cell, errors,
 * note }). Returns handle(message) for the worker's replies, refresh() for when A is ready
 * and reset() for when A is replaced.
 */
export function setupUBRefine(hooks) {
  const dlg = $('ubref');
  // search: the peaks found and the options they were found with; fit: the last refinement.
  let search = null, fit = null, job = 0, busy = false, chosen = null;

  const note = (el, msg, kind = '') => {
    el.textContent = msg;
    el.className = `note ${kind}`;
    el.hidden = !msg;
  };
  const setBusy = (on) => {
    busy = on;
    for (const el of dlg.querySelectorAll('input, select')) el.disabled = on;
    $('ur-run').disabled = $('ur-cancel').disabled = $('ur-close').disabled = on;
  };
  const progress = (label, fraction) => {
    $('ur-progress').hidden = false;
    $('ur-progress-label').textContent = label;
    $('ur-progress-pct').textContent = fraction > 0 ? `${Math.round(100 * fraction)}%` : '';
    $('ur-progress-fill').style.width = `${Math.round(100 * Math.min(1, Math.max(0, fraction)))}%`;
  };

  for (const [key, sys] of Object.entries(UB_SYSTEMS)) $('ur-system').add(new Option(`${sys.label} (${sys.free})`, key));

  /** The panel card's button: available when A has a UB and HKL axes with uniform bins. */
  function refresh() {
    $('ubref-open').disabled = !!reindexProblem(hooks.datasetA());
  }

  function open() {
    const A = hooks.datasetA(), why = reindexProblem(A);
    if (why) return note($('reindex-status'), why, 'error');
    const cell = cellFromUB(A.lattice.ub);
    $('ur-summary').innerHTML = [['Dataset A', esc(hooks.name())], ['Its UB', cellText(cell)]].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    if (!chosen) {
      const laue = hooks.laue(), s = suggestSystem(cell, laue);
      $('ur-system').value = s.system;
      note($('ur-system-note'), s.from === 'laue' ? `From the Laue class ${laue}.` : 'From the metric of A\'s cell (lengths within 1%, angles within 1°).');
    }
    if (!dlg.open) dlg.showModal();
  }
  $('ubref-open').onclick = open;

  // ---- peak search, in A's worker --------------------------------------------------------
  const options = () => ({
    threshold: Number($('ur-threshold').value), maxPeaks: Math.round(Number($('ur-max').value)),
    radius: $('ur-radius').value.trim() ? Math.round(Number($('ur-radius').value)) : null,
  });
  const sameSearch = (o) => search && ['threshold', 'maxPeaks', 'radius'].every((k) => search.options[k] === o[k]);

  $('ur-run').onclick = () => {
    const o = options();
    if (!(o.threshold > 0) || !(o.maxPeaks >= 10) || (o.radius !== null && !(o.radius >= 2 && o.radius <= 20))) {
      return note($('ur-note'), 'Threshold: a positive number; peaks: at least 10; box: 2 to 20 voxels, or empty for automatic.', 'error');
    }
    if (sameSearch(o)) return refit();
    search = fit = null;
    $('ur-result').hidden = true;
    note($('ur-note'), '');
    setBusy(true);
    progress('Finding peaks', 0);
    hooks.post({ type: 'ub-peaks', id: ++job, ...o });
  };

  /** A reply from A's worker; replies to an earlier request are dropped. */
  function handle(msg) {
    if (msg.id !== job) return;
    if (msg.type === 'progress-ub-peaks') {
      progress('Finding peaks', msg.fraction);
    } else if (msg.type === 'ub-peaks') {
      setBusy(false);
      $('ur-progress').hidden = true;
      search = { options: options(), peaks: msg.peaks, stats: msg.stats };
      refit();
    } else if (msg.type === 'ub-peaks-error') {
      setBusy(false);
      $('ur-progress').hidden = true;
      note($('ur-note'), msg.message, 'error');
    }
  }

  // ---- fit, here: fast enough to redo whenever the system or the tolerance changes ----------
  function refit() {
    if (!search) return;
    const A = hooks.datasetA(), ub0 = A.lattice.ub, system = $('ur-system').value, tol = Number($('ur-tol').value);
    if (!(tol > 0 && tol < 0.5)) return note($('ur-note'), 'The index tolerance must be between 0 and 0.5.', 'error');
    const scale = Math.cbrt(voxelQ(A.dims, ub0).reduce((p, x) => p * x, 1));
    try {
      fit = refineUB(search.peaks, ub0, system, { tol, scale });
    } catch (err) {
      fit = null;
      $('ur-result').hidden = true;
      return note($('ur-note'), err.message, 'error');
    }
    let free = null;
    if (system !== 'triclinic') {
      try { free = refineUB(search.peaks, ub0, 'triclinic', { tol, scale }); } catch { /* reported as unavailable */ }
    }
    note($('ur-note'), '');
    showResult(A, ub0, scale, tol, free);
  }
  $('ur-system').onchange = () => { chosen = $('ur-system').value; note($('ur-system-note'), ''); refit(); };
  $('ur-tol').onchange = refit;
  for (const id of ['ur-threshold', 'ur-max', 'ur-radius']) $(id).oninput = () => { $('ur-run').textContent = sameSearch(options()) ? 'Refine' : 'Find peaks and refine'; };

  function showResult(A, ub0, scale, tol, free) {
    const { stats, peaks } = search, r = fit, sys = UB_SYSTEMS[r.system];
    const vox = (q) => fmt(q / scale, 3);
    const box = search.options.radius === null
      ? `± ${stats.radius} voxels, from the peaks' half widths ${stats.halfWidth.map((w) => fmt(w, 2)).join(', ')} voxels` : `± ${stats.radius} voxels`;
    const left = r.within - r.used;
    const shift = indexShift(mul(inv3(r.ub), ub0), A.dims);
    const step = Math.min(...A.dims.map((d) => (d.edges.at(-1) - d.edges[0]) / (d.edges.length - 1)));
    const rows = [
      ['Peaks', `${peaks.length} found (box ${box}); ${r.within} within ${fmt(tol, 3)} of integer indices; ${r.used} used`
        + `${left ? `, ${left} left out with a residual above ${vox(r.cut)} voxel` : ''}`],
      [`Cell (${sys.short})`, `${cellText(r.cell, r.su)}; V ${formatSU(r.volume, r.volumeSU, 2)} Å³`],
      ['Fit', `rms ${fmt(r.rms * 1e3, 3)} × 10⁻³ Å⁻¹ (${vox(r.rms)} voxel); per index ${r.rmsH.map((x) => fmt(x, 4)).join(', ')} r.l.u.; `
        + `${r.rounds} fit${r.rounds === 1 ? '' : 's'}`],
    ];
    if (free) rows.push(['Unconstrained', `${cellText(free.cell, free.su)}; rms ${fmt(free.rms * 1e3, 3)} × 10⁻³ Å⁻¹, ${free.used} peaks`]);
    rows.push(['Orientation', `turned by ${fmt(r.rotation, 4)}° from A's UB`]);
    rows.push(['Reindexing', `moves indices by up to ${fmt(shift, 4)} r.l.u. (${fmt(shift / step, 3)} voxel) at the grid's corners`]);
    const warn = [];
    if (!r.settled) warn.push(`The peaks used changed in each of the ${r.rounds} fits: the result is not stable.`);
    if (r.used < 0.5 * peaks.length) warn.push(`Only ${r.used} of the ${peaks.length} peaks are used.`);
    if (r.rms > 0.5 * scale) warn.push('The residuals average more than half a voxel: the peaks may be split (a volume symmetrized with a wrong UB, twins) or broad.');
    if (free && r.rms > 1.5 * free.rms) warn.push(`The constraint raises the rms by ${fmt(100 * (r.rms / free.rms - 1), 0)} %: the cell may not have this symmetry.`);
    if (warn.length) rows.push(['Check', `<span class="rk-bad">${warn.join(' ')}</span>`]);
    $('ur-report').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('ur-result').hidden = false;
    $('ur-run').textContent = 'Refine';
  }

  // ---- results -------------------------------------------------------------------------------
  const stem = () => hooks.name().replace(/\.[^.]+$/, '');
  const description = () => `refined in the NeXus Viewer from ${fit.used} Bragg peaks of ${hooks.name()}, cell constrained to ${UB_SYSTEMS[fit.system].short}`;
  const errors = () => Object.fromEntries([...KEYS.map((k, i) => [k, fit.su[i]]), ['volume', fit.volumeSU]]);
  const download = (blob, name) => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 10000);
  };
  $('ur-reindex').onclick = () => {
    if (!fit) return;
    const cell = Object.fromEntries(KEYS.map((k, i) => [k, fit.cell[i]]));
    dlg.close();
    hooks.reindexWith(`refined UB (${fit.system}, ${fit.used} peaks)`, { ub: fit.ub, cell, errors: errors(), note: `UB ${description()}` });
  };
  $('ur-save').onclick = () => {
    if (!fit) return;
    download(new Blob([formatIsawUB(fit.ub, { errors: errors(), note: `UB ${description()}.` })], { type: 'text/plain' }), `${stem()}_refined_${fit.system}.mat`);
  };
  $('ur-peaks').onclick = () => {
    if (!fit) return;
    download(new Blob([peaksCSV(search.peaks, fit)], { type: 'text/csv' }), `${stem()}_peaks.csv`);
  };
  $('ur-cancel').onclick = $('ur-close').onclick = () => { if (!busy) dlg.close(); };
  dlg.addEventListener('cancel', (e) => { if (busy) e.preventDefault(); });

  /** Dataset A was replaced: forget its peaks and any reply still on its way. */
  function reset() {
    job++;
    search = fit = chosen = null;
    if (busy) setBusy(false);
    if (dlg.open) dlg.close();
    $('ur-result').hidden = $('ur-progress').hidden = true;
    $('ur-run').textContent = 'Find peaks and refine';
    note($('ur-note'), '');
  }

  return { handle, refresh, reset };
}
