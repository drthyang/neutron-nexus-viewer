import { COLORMAPS } from './colormaps.js';
import { cutAxis, cutBand, cutCrossing, cutEdges } from './cut.js';
import { exportPlan } from './export.js';
import { cartesianBasis, isDirect, nominalCell, planeGeometry, reciprocalMetric } from './nexus.js';
import { parseBins, powderPlan, qExtent } from './powder.js';
import { IDENTITY_MAP } from './slab.js';
import { closeGroup, formatOp, indexMaps, metricChange, parseOps, PRESETS } from './symmetry.js';

const $ = (id) => document.getElementById(id);
// [fixed, x, y] display dimensions: HK, HL and KL planes for Mantid HKL data.
const LAYOUT = [[2, 0, 1], [1, 0, 2], [0, 1, 2]];
const VIEW_KEY = { 2: 'hk', 1: 'hl', 0: 'kl' };
const LUTS = Object.fromEntries(Object.entries(COLORMAPS).map(([name, hex]) =>
  [name, Uint8Array.from(hex.match(/../g), (h) => parseInt(h, 16))]));
const NO_SYMMETRY = { name: '1', ops: [[1, 0, 0, 0, 1, 0, 0, 0, 1]], maps: [IDENTITY_MAP] };
// Canvas colors, matching the page tokens.
const INK = '#121821', INK2 = '#475467', AXIS = '#9aa5b3', MISSING = '#e8ecf1';
const MONO = 'ui-monospace, "SF Mono", Menlo, Consolas, monospace';
const SANS = '-apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, sans-serif';

const svg = (body, size = 15, extra = '') => `<svg width="${size}" height="${size}" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" ${extra}>${body}</svg>`;
const ICONS = {
  save: svg('<path d="M2.5 5.5h2l1.2-2h4.6l1.2 2h2a1 1 0 0 1 1 1v6a1 1 0 0 1-1 1h-11a1 1 0 0 1-1-1v-6a1 1 0 0 1 1-1Z"/><circle cx="8" cy="9.2" r="2.3"/>'),
  focus: svg('<rect x="1.75" y="1.75" width="8" height="12.5" rx="1"/><path d="M11.75 2.5h2.5M11.75 8h2.5M11.75 13.5h2.5"/>'),
  quad: svg('<rect x="1.75" y="1.75" width="5.5" height="5.5" rx="1"/><rect x="8.75" y="1.75" width="5.5" height="5.5" rx="1"/><rect x="1.75" y="8.75" width="5.5" height="5.5" rx="1"/><rect x="8.75" y="8.75" width="5.5" height="5.5" rx="1"/>'),
  max: svg('<path d="M9.5 2.5h4v4M6.5 13.5h-4v-4M13.5 2.5 9 7M2.5 13.5 7 9"/>'),
  restore: svg('<path d="M13.5 6.5h-4v-4M2.5 9.5h4v4M9.5 6.5 14 2M6.5 9.5 2 14"/>'),
  reset: svg('<path d="M2.8 8a5.2 5.2 0 1 0 1.6-3.8"/><path d="M2.5 2.5v3h3"/>'),
  gear: svg('<path d="M2.5 4.5h11M2.5 8h11M2.5 11.5h11"/><circle cx="5.5" cy="4.5" r="1.5" fill="#fff"/><circle cx="10.5" cy="8" r="1.5" fill="#fff"/><circle cx="7" cy="11.5" r="1.5" fill="#fff"/>'),
  download: svg('<path d="M8 2.5v8M4.75 7.25 8 10.5l3.25-3.25"/><path d="M2.5 11v1.5a1 1 0 0 0 1 1h9a1 1 0 0 0 1-1V11"/>'),
  // How compared datasets share a slice: split along the diagonal (two), or quadrants (three or four).
  diagonal: svg('<path d="M1.75 1.75v12.5h12.5Z" fill="currentColor" opacity="0.3" stroke="none"/><rect x="1.75" y="1.75" width="12.5" height="12.5" rx="1"/><path d="M2.25 2.25l11.5 11.5"/>', 14),
  quadrants: svg('<path d="M1.75 8h6.25v6.25H1.75Z" fill="currentColor" opacity="0.3" stroke="none"/><rect x="1.75" y="1.75" width="12.5" height="12.5" rx="1"/><path d="M8 1.75v12.5M1.75 8h12.5"/>', 14),
};

let worker = null, meta = null, panels = [], settings = null, sourceName = '', sourceSize = 0, autoscaled = false;
// The colormap last chosen (and saved); signed data opens with coolwarm instead.
let cmapChoice = 'magma';
let requestId = 0, symmetry = NO_SYMMETRY, mask = null;
// Workspace layout: 'quad', 'focus' or 'single', around the primary view ('hk', 'hl', 'kl' or '3d').
let layout = 'quad', primary = 'hk', lastMulti = 'quad';
const views = {};
// What clicking a slice does: 'navigate' moves the other two slices, 'zoom' zooms in,
// 'move' drags the visible region and 'cut' draws a line cut.
let clickMode = 'navigate';
// 3-D view state: the lazily loaded View3D, its elements and the isosurface request queue.
let view3d = null, iso = null;
// I(Q) and line cut views: their elements, zoom and one request queue per dataset
// (`a` for A, `more[i]` for the dataset in slot i).
// The 3-D view, I(Q) and the cut share the fourth place in the layouts; `slot` is the one shown.
let powder = null, cut = null, slot = '3d', powderScale = 'linear', cutScale = 'linear';
const SLOT_VIEWS = ['3d', 'iq', 'cut'];
// Counts mask rebuilds, so I(Q) is recomputed after each.
let maskVersion = 0;
// Datasets compared with A: B, C and D in slots 0, 1 and 2 (null when empty), each
// with its own worker, index maps and mask. `compareView` is what the slices show:
// 'split' (all of them: two split along the diagonal, A below; three or four in
// quadrants, A, B, C, D from lower left to upper right), or one alone, 'a' … 'd'.
// `pendingCompare` holds the URLs to open as B, C and D once A is ready.
let slots = [null, null, null], compareView = 'split', pendingCompare = [];
const LETTERS = ['A', 'B', 'C', 'D'];
// Symmetry and mask to apply once the next file opens (?sym= and ?mask=, or the demo).
let pendingProcessing = null;
const DEMO = { url: 'examples/demo_300K.nxs', compare: 'examples/demo_10K.nxs', sym: '6/mmm', mask: '1' };
// True while the color range is the automatic one (not edited by hand).
let rangeIsAuto = false;
// The export for NEBULA3D being built: { name, send, attrs }. `handoff` is the
// NEBULA3D tab it is sent to (Open in NEBULA3D): { id, win | channel, origin, ready, file, sent }.
let exportJob = null, handoff = null;
// ?nebula3d= points Open in NEBULA3D at another deployment (a local dev server).
const NEBULA3D_URL = new URLSearchParams(location.search).get('nebula3d') || 'https://drthyang.github.io/nebula3d/';
const newLayer = () => ({ data: null, version: 0, busy: false, wanted: null, image: null, imageKey: null, error: '' });
const others = () => slots.filter(Boolean);
const readyOthers = () => slots.filter((d) => d?.ready);
/** What the slices show while comparing ('split', or the letter of one dataset that is ready), or null. */
function compareShown() {
  if (!readyOthers().length) return null;
  const k = 'abcd'.indexOf(compareView);
  return k > 0 && !slots[k - 1]?.ready ? 'split' : compareView;
}
/** The dataset shown alone (0 for A, 1–3 for B–D), or 0 while showing several. */
const shownIndex = () => Math.max(0, 'abcd'.indexOf(compareShown()));
/** Dataset k (0 for A) of a panel or plot: A's layer is the panel itself (a plot's `a`), the others in `more`. */
const layerOf = (obj, k) => (k ? obj.more[k - 1] : obj.a ?? obj);
const workerOf = (k) => (k ? slots[k - 1].worker : worker);
/** Letters or names in running text: "A", "A and B", "A, B and C". */
const listText = (items) => (items.length > 1 ? `${items.slice(0, -1).join(', ')} and ${items.at(-1)}` : items[0] ?? '');

// ---- Formatting -------------------------------------------------------------

const sig = (x, n = 2) => Number(x.toPrecision(n));
const roundTo = (x, step) => Number((Math.round(x / step) * step).toPrecision(12));
const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));
const fmt = (x, digits = 3) => String(Number(x.toFixed(digits))).replace('-', '−');
const mb = (bytes) => `${(bytes / 1e6).toFixed(bytes < 1e8 ? 1 : 0)} MB`;
const pct = (x) => `${(100 * x).toFixed(x < 0.01 ? 2 : 1)}%`;
const withUnits = (dim) => (dim.units ? `${dim.label} (${dim.units})` : dim.label);
const escapeHTML = (s) => String(s).replace(/[&<>"']/g, (ch) => `&#${ch.charCodeAt(0)};`);

function fmtValue(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  const s = a >= 1e5 || a < 1e-3 ? v.toExponential(1).replace(/\.0e/, 'e').replace('e+', 'e') : String(sig(v, 3));
  return s.replace('-', '−');
}

function niceStep(raw) {
  const p = 10 ** Math.floor(Math.log10(raw));
  return [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw * 0.999);
}

function niceTicks(lo, hi, count = 5, fixedStep = null) {
  const raw = (hi - lo) / count;
  if (!(raw > 0)) return [];
  const step = fixedStep ?? niceStep(raw);
  const out = [];
  for (let k = Math.ceil(lo / step - 1e-9); k * step <= hi + step * 1e-9; k++) out.push(roundTo(k * step, step / 100));
  return out;
}

// ---- Small widgets -------------------------------------------------------------

/** Paint the filled part of a range input (nebula3d style, via --p). */
function paint(range) {
  const lo = Number(range.min || 0), hi = Number(range.max || 100);
  range.style.setProperty('--p', `${hi > lo ? (100 * (Number(range.value) - lo)) / (hi - lo) : 0}%`);
}
const paintAll = () => document.querySelectorAll('input[type=range]').forEach(paint);
document.addEventListener('input', (e) => { if (e.target.type === 'range') paint(e.target); });

/** Segmented control: value in data-value; fires onChange after a click. */
function segmented(el, onChange) {
  el.addEventListener('click', (e) => {
    const button = e.target.closest('button');
    if (!button || button.classList.contains('on')) return;
    setSegmented(el, button.dataset.value);
    onChange(button.dataset.value);
  });
}
function setSegmented(el, value) {
  el.dataset.value = value;
  for (const b of el.querySelectorAll('button')) b.classList.toggle('on', b.dataset.value === value);
}

/**
 * `name` shortened with "…" until `fits` accepts it. On its own, a name keeps
 * its start and its last 12 characters (often a temperature or run number).
 * Next to the name of the file it is compared with (`other`), it keeps the part
 * where the two names differ, with as much context around it as fits.
 */
function shortenName(name, fits, other = null) {
  if (fits(name)) return name;
  const label = (s, e) => `${s > 0 ? '…' : ''}${name.slice(s, e)}${e < name.length ? '…' : ''}`;
  if (!other || other === name) {
    const tail = name.slice(-12);
    let head = name.length - tail.length;
    while (head > 0 && !fits(`${name.slice(0, head)}…${tail}`)) head--;
    return `${name.slice(0, head)}…${tail}`;
  }
  let p = 0, q = 0;
  while (p < name.length && name[p] === other[p]) p++;
  while (q < name.length - p && q < other.length - p && name.at(-1 - q) === other.at(-1 - q)) q++;
  let s = Math.min(p, name.length - 1), e = Math.max(s + 1, name.length - q);
  while (e > s + 1 && !fits(label(s, e))) e--;
  // Widen the kept part, toward the start first, while it fits.
  for (let grew = true; grew;) {
    grew = false;
    if (s > 0 && fits(label(s - 1, e))) { s--; grew = true; }
    if (e < name.length && fits(label(s, e + 1))) { e++; grew = true; }
  }
  return label(s, e);
}

/** A file name for running text: at most `max` characters, shortened in the middle. */
const shortName = (name, max = 36) => shortenName(name, (t) => t.length <= max);

// File names in the page are fitted to their element (CSS .fname), again whenever it resizes.
const nameFitter = new ResizeObserver((entries) => { for (const { target } of entries) fitName(target); });

/** Show a file name in `el`, shortened to fit (see shortenName); the full name is the tooltip. */
function setFileName(el, name, other = null) {
  el.dataset.name = name;
  if (other) el.dataset.other = other;
  else delete el.dataset.other;
  el.title = name;
  fitName(el);
  nameFitter.observe(el);
}

function fitName(el) {
  const fits = (text) => {
    el.textContent = text;
    return el.scrollWidth <= el.clientWidth;
  };
  el.textContent = shortenName(el.dataset.name ?? '', fits, el.dataset.other ?? null);
}

const urlName = (url) => decodeURIComponent(new URL(url, location.href).pathname.split('/').pop()) || 'remote.nxs';

// Work in progress, shown in the top bar: building a mask, I(Q), the NEBULA3D
// export, opening B. The status names the job that started first and how far
// it is (with the number of others), the bar under the top bar shows all of
// them together, and hovering lists them. Without jobs, the status is the one
// last set with status().
const jobs = new Map();
let idle = ['', 'no file'];

function status(state, text) {
  idle = [state, text];
  showJobs();
}

/** Start or update job `key`: `name` for the top bar, `step` for its tooltip. */
function setJob(key, name, step = '', fraction = 0) {
  jobs.set(key, { name, step, fraction: clamp(fraction || 0, 0, 1) });
  showJobs();
}

function endJob(...keys) {
  for (const key of keys) jobs.delete(key);
  showJobs();
}

function showJobs() {
  const list = [...jobs.values()], busy = list.length > 0, first = list[0];
  const percent = (f) => `${Math.round(100 * f)}%`;
  $('status-dot').className = `status-dot ${busy ? 'busy' : idle[0]}`;
  $('status-text').textContent = busy ? `${first.name} · ${percent(first.fraction)}${list.length > 1 ? ` +${list.length - 1}` : ''}` : idle[1];
  $('status').title = busy
    ? list.map((j) => `${j.name}${j.step ? `: ${j.step}` : ''} · ${percent(j.fraction)}`).join('\n')
    : 'Everything runs in this browser tab';
  $('work-progress').hidden = !busy;
  if (busy) $('work-progress-fill').style.width = percent(list.reduce((s, j) => s + j.fraction, 0) / list.length);
}

// Popovers: one open at a time, placed under the button that opened it.
function closePopovers() {
  for (const pop of document.querySelectorAll('.popover')) pop.hidden = true;
  for (const b of document.querySelectorAll('[data-pop]')) b.setAttribute('aria-expanded', 'false');
}

function togglePopover(button) {
  const pop = $(button.dataset.pop), open = pop.hidden;
  closePopovers();
  if (!open) return;
  pop.hidden = false;
  button.setAttribute('aria-expanded', 'true');
  const r = button.getBoundingClientRect(), w = pop.offsetWidth;
  const left = r.left + r.width / 2 < innerWidth / 2 ? r.left : r.right - w;
  pop.style.left = `${clamp(left, 12, innerWidth - w - 12)}px`;
  pop.style.top = `${r.bottom + 8}px`;
}

document.addEventListener('click', (e) => {
  const button = e.target.closest('[data-pop]');
  if (button) togglePopover(button);
  else if (!e.target.closest('.popover')) closePopovers();
});
addEventListener('resize', closePopovers);
addEventListener('keydown', (e) => {
  if (e.key !== 'Escape') return;
  if (document.querySelector('.popover:not([hidden])')) closePopovers();
  else if (layout === 'single') setLayout(lastMulti, primary);
});

// ---- Files and the worker ---------------------------------------------------

function show(stage) {
  $('intro').hidden = stage !== 'intro';
  $('loading').hidden = stage !== 'loading';
  $('app').hidden = stage !== 'workspace';
  for (const el of document.querySelectorAll('[data-loaded]')) el.hidden = stage !== 'workspace';
}

/** Loading card: the current step, and the fraction done (0 hides the percentage). */
function progress(label, fraction) {
  $('progress-label').textContent = label;
  $('progress-pct').textContent = fraction > 0 ? `${Math.round(100 * clamp(fraction, 0, 1))}%` : '';
  $('progress').style.width = `${Math.round(100 * clamp(fraction, 0, 1))}%`;
}

function error(message) {
  $('error-text').textContent = message;
  $('error').hidden = !message;
}

function fail(message) {
  error(message);
  endJob('mask-a', 'iq-a', 'export');
  status(meta ? 'ok' : '', meta ? 'ready' : 'no file');
  if (!panels.length) show('intro');
}

function openFile(file) {
  worker?.terminate();
  for (const d of others()) d.worker?.terminate();
  slots = [null, null, null];
  exportJob = null;
  jobs.clear();
  showCompare();
  closePopovers();
  meta = null;
  panels = [];
  autoscaled = false;
  symmetry = NO_SYMMETRY;
  mask = null;
  iso = null;
  powder = null;
  sourceName = file.name;
  sourceSize = file.size;
  $('workspace').replaceChildren();
  error('');
  show('loading');
  status('busy', 'reading');
  setFileName($('dataset-name'), file.name);
  $('dataset-facts').hidden = true;
  $('loading-title').textContent = 'Opening file';
  setFileName($('loading-name'), file.name);
  $('loading-size').textContent = mb(file.size);
  progress('Starting', 0);
  document.title = `${file.name} · NeXus Viewer`;
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = ({ data }) => handlers[data.type]?.(data);
  worker.onerror = (e) => {
    e.preventDefault();
    fail(e.message || 'The HDF5 reader could not start. It is loaded from cdn.jsdelivr.net; check the network connection.');
  };
  worker.postMessage({ type: 'open', file });
}

/** Download `url` into a File, reporting (bytes so far, total or 0). */
async function fetchFile(url, onProgress) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const total = Number(res.headers.get('content-length')) || 0;
  const reader = res.body.getReader(), parts = [];
  for (let got = 0; ;) {
    const { done, value } = await reader.read();
    if (done) break;
    parts.push(value);
    got += value.length;
    onProgress(got, total);
  }
  return new File(parts, urlName(url));
}

const downloaded = (got, total) => `${mb(got)}${total ? ` of ${mb(total)}` : ''}`;

async function openURL(url) {
  show('loading');
  status('busy', 'downloading');
  $('loading-title').textContent = 'Downloading';
  setFileName($('loading-name'), urlName(url));
  $('loading-size').textContent = '';
  progress(`from ${new URL(url, location.href).host}`, 0);
  try {
    openFile(await fetchFile(url, (got, total) => progress(downloaded(got, total), total ? got / total : 0)));
  } catch (err) {
    fail(`Could not download ${url}: ${err.message}. The server must allow cross-origin requests.`);
  }
}

const handlers = {
  progress: ({ label, fraction }) => progress(label, fraction),
  meta: ({ info }) => { meta = info; },
  ready: ({ stats, seconds }) => {
    meta.stats = stats;
    meta.seconds = seconds;
    status('ok', 'ready');
    setupViewer();
    show('workspace');
    showCompare();
    redraw();
    panels.forEach(request);
    setup3D();
    for (const url of pendingCompare) openCompareURL(url);
    pendingCompare = [];
    applyPendingProcessing();
  },
  slice: (msg) => {
    const p = panels.find((q) => q.fixed === msg.fixed);
    if (!p) return;
    p.busy = false;
    p.data = msg;
    p.version++;
    showCaption(p);
    if (p.wanted) send(p);
    if (!autoscaled && panels.every((q) => q.data)) {
      autoscaled = true;
      autoRange();
    }
    redraw();
  },
  error: ({ fixed, message }) => {
    const p = panels.find((q) => q.fixed === fixed);
    if (!p) return fail(message);
    p.busy = false;
    p.caption.textContent = message;
    p.caption.title = message;
    if (p.wanted) send(p);
  },
  iso: (msg) => {
    iso.busy = false;
    iso.level = msg.level;
    iso.range = msg.range;
    if (!iso.userLevel) iso.levelInput.value = sig(msg.level, 3);
    iso.slider.value = levelToSlider(msg.level);
    paint(iso.slider);
    const triangles = msg.indices.length / 3;
    iso.caption.textContent = `${msg.shape[0]}³ grid · ${triangles.toLocaleString()} tri${msg.symmetry !== '1' ? ` · ${msg.symmetry}` : ''}${mask ? ' · masked' : ''}`;
    iso.caption.title = `${msg.shape.join(' × ')} blocks (${msg.factor}× binned), level ${fmtValue(msg.level)}, ${triangles.toLocaleString()} triangles${msg.note ? `. ${msg.note}` : ''}`;
    view3d.setMesh(msg.positions, msg.indices, Number(iso.opacity.value));
    if (iso.wanted) sendIso();
  },
  'iso-error': ({ message }) => {
    iso.busy = false;
    iso.caption.textContent = message;
    if (iso.wanted) sendIso();
  },
  'progress-mask': ({ label, fraction }) => {
    setJob('mask-a', 'Mask', label, fraction);
    $('mask-status').className = 'note';
    $('mask-status').textContent = `${label}… ${Math.round(100 * fraction)}%`;
  },
  mask: ({ stats, radius, k, symmetry: group, seconds }) => {
    endJob('mask-a');
    $('mask-apply').disabled = false;
    mask = stats ? { ...stats, radius, k, group } : null;
    maskVersion++;
    // Other datasets may have opened while this mask was being built.
    for (const d of readyOthers()) if (d.maskRequested !== `${radius},${k}`) sendCompareMask(d, radius, k);
    $('mask-clear').disabled = $('mask-download').disabled = $('mask-removed').disabled = !mask;
    if (!mask) $('mask-removed').checked = false;
    showMask(seconds);
    panels.forEach((p) => request(p, 'a'));
    requestIso();
    requestPowder();
    requestCut();
    describe();
  },
  'mask-error': ({ message }) => {
    endJob('mask-a');
    $('mask-apply').disabled = false;
    $('mask-status').className = 'note error';
    $('mask-status').textContent = message;
    updateStates();
  },
  'mask-file': ({ blob }) => download(blob, `${stem()}_mask.npy.gz`),
  'progress-export': ({ label, fraction }) => {
    showExportProgress(label, fraction);
    forwardProgress(label, fraction);
  },
  'export-file': ({ blob, stats, seconds }) => {
    const { name, send, attrs } = exportJob;
    exportJob = null;
    endJob('export');
    showExport();
    if (send) {
      if (!handoff) return; // the NEBULA3D tab was closed meanwhile
      handoff.file = new File([blob], name, { type: 'application/x-hdf5' });
      handoff.meta = attrs;
      exportNote('busy', 'sending…', `Built ${shortName(name)} (${mb(blob.size)}); waiting for NEBULA3D to start…`, name);
      sendHandoff();
      return;
    }
    download(blob, name);
    exportNote('ok', 'saved', `Saved ${shortName(name)} (${mb(blob.size)}, ${pct(stats.valid / stats.total)} of voxels valid) in ${seconds.toFixed(1)} s. `
      + 'In NEBULA3D, open it with Load volume…', name);
  },
  'export-error': ({ message }) => {
    if (exportJob?.send) endHandoff(`the viewer could not build the volume: ${message}`);
    exportJob = null;
    endJob('export');
    showExport();
    $('export-status').className = 'note error';
    $('export-status').textContent = message;
  },
  'progress-powder': (msg) => powderProgress(0, msg),
  powder: (msg) => powderResult(0, msg),
  'powder-error': ({ message }) => powderResult(0, null, message),
  'progress-cut': (msg) => cutProgress(0, msg),
  cut: (msg) => cutResult(0, msg),
  'cut-error': ({ message }) => cutResult(0, null, message),
};

/**
 * Ask for the panel's slice of dataset 'a', of one other dataset (its slot
 * object), or of 'both' (every dataset that is loaded) with the cut.
 */
function request(p, which = 'both') {
  const center = Number(p.center.value), thickness = Number(p.width.value);
  if (!Number.isFinite(center) || !Number.isFinite(thickness) || thickness <= 0) {
    p.caption.textContent = 'Enter a finite center and a positive thickness.';
    return;
  }
  const one = typeof which === 'object' ? which : null;
  if (!one) {
    p.wanted = { center, thickness };
    if (!p.busy) send(p);
  }
  for (const d of which === 'a' ? [] : readyOthers()) {
    if (one && one !== d) continue;
    const layer = p.more[d.slot];
    layer.wanted = { center, thickness };
    if (!layer.busy) sendOther(p, d);
  }
  // The cut lies in its slice's plane and follows it.
  if (which === 'both' && cutPanel() === p) requestCut();
}

// One request in flight per panel; slider drags coalesce to the latest value.
function send(p) {
  const q = p.wanted;
  p.wanted = null;
  p.busy = true;
  worker.postMessage({
    type: 'slice', id: ++requestId, fixed: p.fixed, ...q,
    maps: symmetry.maps, symmetry: symmetry.name, removed: $('mask-removed').checked,
  });
}

function sendOther(p, d) {
  const layer = p.more[d.slot], q = layer.wanted;
  layer.wanted = null;
  layer.busy = true;
  d.worker.postMessage({
    type: 'slice', id: ++requestId, fixed: p.fixed, ...q,
    maps: d.maps, symmetry: d.maps.length > 1 ? symmetry.name : '1', removed: $('mask-removed').checked,
  });
}

/** View header: position, bins and coverage, and the other datasets' coverage when comparing. */
function showCaption(p) {
  const F = meta.dims[p.fixed], d = p.data;
  if (!d) return;
  p.pos.textContent = `${F.label} = ${fmt(d.center)}`;
  const slab = (x) => `${F.label} ∈ [${fmt(x.slab[0])}, ${fmt(x.slab[1])}], ${x.bins} bin${x.bins === 1 ? '' : 's'}, ${pct(x.coverage)} of the plane measured`
    + `${x.order > 1 ? `, averaged over ${x.symmetry} (${x.order} operations)` : ''}${x.removed ? ', removed voxels only' : ''}`;
  let text = `${d.bins} bin${d.bins === 1 ? '' : 's'} · ${pct(d.coverage)}${d.order > 1 ? ` · ${d.symmetry}` : ''}${d.removed ? ' · removed only' : ''}`;
  let title = slab(d);
  const more = readyOthers();
  if (more.length) {
    const layers = more.map((o) => [o.letter, p.more[o.slot]]);
    text = `A ${text} | ${layers.map(([k, l]) => `${k} ${l.data ? pct(l.data.coverage) : l.error ? 'no data' : '…'}`).join(' | ')}`;
    title = `A: ${title}\n${layers.map(([k, l]) => `${k}: ${l.data ? slab(l.data) : l.error || 'loading'}`).join('\n')}`;
  }
  p.caption.textContent = text;
  p.caption.title = title;
}

const symmetryNote = (data) => (data.order > 1 ? ` · ${data.symmetry} (${data.order} ops)` : '');
const stemOf = (name) => name.replace(/\.[^.]+$/, '');
const stem = () => stemOf(sourceName);

/**
 * File stems as one name: the words they all share once, and the ones that
 * differ as "A_vs_B_vs_C" (demo_300K and demo_10K give demo_300K_vs_10K).
 */
function setStem(stems) {
  if (stems.length < 2) return stems[0] ?? '';
  const words = (s) => s.match(/[^_\-. ]+|[_\-. ]+/g) ?? [];
  const trim = (parts) => parts.join('').replace(/^[_\-. ]+|[_\-. ]+$/g, '');
  const ws = stems.map(words), n = Math.min(...ws.map((w) => w.length));
  let head = 0, tail = 0;
  while (head < n && ws.every((w) => w[head] === ws[0][head])) head++;
  while (tail < n - head && ws.every((w) => w[w.length - 1 - tail] === ws[0][ws[0].length - 1 - tail])) tail++;
  const mids = ws.map((w) => trim(w.slice(head, w.length - tail)));
  if (mids.some((m) => !m)) return stems.join('_vs_');
  return [trim(ws[0].slice(0, head)), mids.join('_vs_'), trim(ws[0].slice(ws[0].length - tail))].filter(Boolean).join('_');
}

// ---- Viewer setup -------------------------------------------------------------

function describe() {
  const { dims, shape, lattice, stats } = meta;
  const bins = [0, 1, 2].map((d) => shape[2 - d]);
  const widths = dims.map((d) => fmt((d.edges[d.edges.length - 1] - d.edges[0]) / (d.edges.length - 1), 4));

  // Dataset A in the top bar: one line of facts, the full ones as its tooltip.
  const { a, b, c, alpha, beta, gamma } = lattice ?? {};
  const recip = recipOf(lattice);
  $('dataset-facts').textContent = datasetFacts(meta);
  $('dataset-facts').hidden = false;
  $('open').title = `${others().length ? 'Dataset A: ' : ''}${sourceName} (${mb(sourceSize)})\n${datasetDetails(meta, mask)}\nClick to open another file.`;
  $('pipe-measured').textContent = `${(stats.valid / 1e6).toFixed(1)} M · ${pct(stats.fraction)} of the grid`;

  // Full details in the info popover.
  const items = [
    ['File', `${sourceName} (${mb(sourceSize)})`],
    ['Signal', meta.signal],
    ['Errors', meta.errors ? `${meta.errors.path}${meta.errors.squared ? '' : ' (σ, squared when read)'}` : 'none in the file'],
    ['Axes', dims.map((d) => `${d.longName} ${fmt(d.edges[0], 2)} … ${fmt(d.edges[d.edges.length - 1], 2)}`).join('\n')],
    ['Grid', `${bins.join(' × ')} bins, Δ = ${widths.join(', ')}`],
    ['Measured', `${pct(stats.fraction)} of voxels · read in ${meta.seconds.toFixed(1)} s`],
  ];
  if (lattice) {
    items.push([`Cell (${lattice.source})`, `a ${a.toFixed(4)}, b ${b.toFixed(4)}, c ${c.toFixed(4)} Å\nα ${alpha.toFixed(3)}°, β ${beta.toFixed(3)}°, γ ${gamma.toFixed(3)}°`]);
    const [as, bs, cs] = recip.len, [al, be, ga] = recip.angles;
    items.push(['Reciprocal', `a* ${as.toFixed(5)}, b* ${bs.toFixed(5)}, c* ${cs.toFixed(5)} Å⁻¹ (no 2π)\nα* ${al.toFixed(3)}°, β* ${be.toFixed(3)}°, γ* ${ga.toFixed(3)}°`]);
    if ($('angles').checked) items.push(['Display', 'Nominal angles: direct-cell angles within 1° of 60/90/120° are snapped for drawing']);
  }
  items.push(['Symmetry', symmetry.ops.length > 1 ? `${symmetry.name} (${symmetry.ops.length} operations)` : 'none']);
  items.push(['Mask', mask ? `${pct((mask.edge + mask.outlier) / mask.measured)} of voxels removed` : 'none']);
  for (const o of readyOthers()) {
    const M = o.meta, l = M.lattice;
    items.push([`Compare (${o.letter})`, `${o.name} (${mb(o.size)})\n${[0, 1, 2].map((d) => M.shape[2 - d]).join(' × ')} bins, ${pct(M.stats.fraction)} measured`
      + (l ? `\na ${l.a.toFixed(4)}, b ${l.b.toFixed(4)}, c ${l.c.toFixed(4)} Å, γ ${l.gamma.toFixed(3)}°` : '')
      + (o.mask ? `\nmask: ${pct((o.mask.edge + o.mask.outlier) / o.mask.measured)} removed` : '')]);
  }
  $('info-meta').innerHTML = items.map(([k, v]) => `<dt>${escapeHTML(k)}</dt><dd>${escapeHTML(v).replace(/\n/g, '<br>')}</dd>`).join('');
  updateStates();
  showExport();
}

/** Reciprocal lengths (Å⁻¹, no 2π) and angles (α*, β*, γ*) of a cell, or null. */
function recipOf(lattice) {
  if (!lattice) return null;
  const G = reciprocalMetric(lattice), len = [0, 1, 2].map((i) => Math.sqrt(G[i][i]));
  const ang = (i, j) => Math.acos(G[i][j] / (len[i] * len[j])) * 180 / Math.PI;
  return { len, angles: [ang(1, 2), ang(0, 2), ang(0, 1)] };
}

/** A dataset's cell (angles other than 90°), grid and measured fraction, in one short line for the top bar. */
function datasetFacts(m) {
  const l = m.lattice, bins = [0, 1, 2].map((d) => m.shape[2 - d]);
  const same = (xs) => xs.every((x) => x === xs[0]);
  const parts = [];
  if (l) {
    // Angles only where they differ from 90°.
    const angles = [['α', l.alpha], ['β', l.beta], ['γ', l.gamma]].filter(([, x]) => Math.abs(x - 90) >= 0.05);
    parts.push(`${[l.a, l.b, l.c].map((x) => x.toFixed(3)).join(' ')} Å${angles.map(([k, x], i) => `${i ? '' : ','} ${k} ${Number(x.toFixed(1))}°`).join('')}`);
  }
  parts.push(same(bins) ? `${bins[0]}³` : bins.join('×'), `${pct(m.stats.fraction)} measured`);
  return parts.join(' · ');
}

/** A dataset's cell, reciprocal lattice, grid, measured fraction and mask, one per line (for tooltips). */
function datasetDetails(m, userMask) {
  const l = m.lattice, r = recipOf(l), bins = [0, 1, 2].map((d) => m.shape[2 - d]);
  const widths = m.dims.map((d) => fmt((d.edges[d.edges.length - 1] - d.edges[0]) / (d.edges.length - 1), 4));
  return [
    l ? `Cell ${l.a.toFixed(4)} ${l.b.toFixed(4)} ${l.c.toFixed(4)} Å, ${l.alpha.toFixed(2)}° ${l.beta.toFixed(2)}° ${l.gamma.toFixed(2)}° (${l.source})` : 'No unit cell in the file',
    r ? `Reciprocal ${r.len.map((x) => x.toFixed(4)).join(' ')} Å⁻¹ (no 2π)` : null,
    `Grid ${bins.join(' × ')}, Δ ${widths.join(' ')}`,
    `Measured ${pct(m.stats.fraction)} of voxels${userMask ? `, ${pct((userMask.edge + userMask.outlier) / userMask.measured)} masked` : ''}`,
  ].filter(Boolean).join('\n');
}

/** Processing pipeline and section badges. */
function updateStates() {
  const sym = symmetry.ops.length > 1;
  $('sym-state').textContent = sym ? `${symmetry.name} · ${symmetry.ops.length}` : 'none';
  $('sym-state').className = `state${sym ? ' on' : ''}`;
  $('symmetry').classList.toggle('active', sym);
  $('pipe-sym').className = sym ? 'on' : 'off';
  $('pipe-sym-text').textContent = sym ? `${symmetry.name} · ${symmetry.ops.length} operations` : 'none — voxels used as measured';
  const removed = mask ? pct((mask.edge + mask.outlier) / mask.measured) : null;
  if (!$('mask-apply').disabled) {
    $('mask-state').textContent = mask ? `${removed} removed` : 'off';
    $('mask-state').className = `state${mask ? ' ok' : ''}`;
  }
  $('mask').classList.toggle('active', !!mask);
  $('mask').classList.toggle('mask-on', !!mask);
  $('pipe-mask').className = mask ? 'ok' : 'off';
  $('pipe-mask-text').textContent = mask
    ? `${removed} removed${mask.radius ? ` · edge ${mask.radius}` : ''}${mask.k ? ` · ${mask.k}σ outliers` : ''}${$('mask-removed').checked ? ' · showing removed' : ''}`
    : 'off — all measured voxels';
  // One-line summaries shown on collapsed panel sections.
  const shown = ['A', ...readyOthers().map((o) => o.letter)];
  $('pipe-views').textContent = shown.length > 1
    ? `3 slices, ${shown.join(' | ')} ${others().length > 1 ? 'quadrants' : 'split'} + 3-D (A) + I(Q)` : '3 slices + 3-D + I(Q)';
  $('sum-processing').textContent = `${mask ? `mask ${removed}` : 'no mask'} · ${sym ? symmetry.name : 'no symmetry'}`;
}

function viewShell(key, badge, title) {
  const section = document.createElement('section');
  section.className = 'view';
  section.dataset.key = key;
  section.innerHTML = `
    <header class="view-head">
      ${badge}<span class="view-title">${title}</span><span class="view-pos"></span>
      <span class="view-caption"></span>
      <span class="view-actions"></span>
    </header>
    <div class="view-body"></div>`;
  $('workspace').append(section);
  const actions = section.querySelector('.view-actions');
  const focusBtn = document.createElement('button'), maxBtn = document.createElement('button');
  focusBtn.className = maxBtn.className = 'icon-btn layout-btn';
  focusBtn.type = maxBtn.type = 'button';
  focusBtn.onclick = () => setLayout(layout === 'focus' && primary === key ? 'quad' : 'focus', key);
  maxBtn.onclick = () => setLayout(layout === 'single' && primary === key ? lastMulti : 'single', key);
  section.querySelector('.view-head').ondblclick = (e) => { if (!e.target.closest('button')) maxBtn.click(); };
  section.querySelector('.slot-switch')?.addEventListener('click', (e) => {
    const button = e.target.closest('button');
    if (!button || button.classList.contains('on')) return;
    if (button.dataset.value === 'iq') powder.requested = true;
    setSlot(button.dataset.value);
  });
  // In the focus layout the small views are thumbnails: a click (outside their buttons) shows one large.
  const enlarge = () => section.classList.contains('thumb');
  section.addEventListener('click', (e) => { if (enlarge() && !e.target.closest('button, input, select, a')) setLayout('focus', key); });
  section.addEventListener('keydown', (e) => {
    if (e.target === section && (e.key === 'Enter' || e.key === ' ') && enlarge()) {
      e.preventDefault();
      setLayout('focus', key);
    }
  });
  views[key] = { section, focusBtn, maxBtn };
  return { section, actions, focusBtn, maxBtn, q: (sel) => section.querySelector(sel) };
}

/** The badge of the 3-D view, I(Q) and the cut, which share the fourth place: a switch between them, `key` on. */
const slotSwitch = (key) => `<div class="segmented slot-switch" role="group" aria-label="Fourth view">${
  [['3d', '3D', 'Show the 3-D isosurface here'], ['iq', 'I(Q)', 'Show I(Q), the powder average, here'], ['cut', 'Cut', 'Show the line cut here']].map(([k, label, title]) =>
    `<button type="button" data-value="${k}"${k === key ? ' class="on"' : ` title="${title}"`}>${label}</button>`).join('')}</div>`;

function setupViewer() {
  // A diverging colormap for signed data (a ΔPDF), the chosen one otherwise.
  $('cmap').value = meta.signed ? 'coolwarm' : cmapChoice;
  describe();
  $('angles-wrap').hidden = !LAYOUT.some(([, x, y]) => planeGeometry(meta.dims, meta.lattice, x, y).lattice);
  // Placeholders until the first three slices arrive and autoRange() runs.
  $('vmin').value = 0;
  $('vmax').value = 1;
  $('soft').value = 0.05;
  $('sym-preset').value = '1';
  $('sym-ops').value = '';
  showSymmetry();
  $('mask-clear').disabled = $('mask-download').disabled = $('mask-removed').disabled = true;
  $('mask-removed').checked = false;
  showMask();

  for (const [fixed, x, y] of LAYOUT) {
    const F = meta.dims[fixed], X = meta.dims[x], Y = meta.dims[y], e = F.edges, n = e.length - 1;
    const step = sig((e[n] - e[0]) / n, 2);
    const lo = (e[0] + e[1]) / 2, hi = (e[n - 1] + e[n]) / 2;
    const center = lo <= 0 && hi >= 0 ? 0 : roundTo((lo + hi) / 2, step);
    const maxWidth = Math.min(e[n] - e[0], 41 * step);
    const key = VIEW_KEY[fixed];
    const shell = viewShell(key, `<span class="badge hue-${fixed}">${escapeHTML(F.label)}</span>`, `${escapeHTML(X.label)} – ${escapeHTML(Y.label)}`);
    shell.q('.view-pos').insertAdjacentHTML('beforebegin', '<span class="view-angle" hidden></span>');
    shell.actions.innerHTML = `
      <button type="button" class="btn btn-ghost btn-xs zoom-reset" hidden title="Back to the full view (or double-click the plot)">Reset zoom</button>
      <button type="button" class="icon-btn save" title="Save PNG">${ICONS.save}</button>`;
    shell.actions.append(shell.focusBtn, shell.maxBtn);
    shell.q('.view-body').innerHTML = `<canvas class="plot" role="img"></canvas><span class="overlay-chip" hidden></span>`;
    const units = F.units ? escapeHTML(F.units) : '';
    const foot = document.createElement('footer');
    foot.className = `view-foot cut--${fixed}`;
    foot.innerHTML = `
      <span class="foot-label" title="Slab center">${escapeHTML(F.label)}</span>
      <input class="slider foot-grow" type="range" aria-label="${escapeHTML(F.label)} center">
      <input class="center num" type="number" step="any" aria-label="${escapeHTML(F.label)} center">
      <span class="unit">${units}</span>
      <span class="foot-sep"></span>
      <span class="foot-label" title="Full slab thickness">Δ${escapeHTML(F.label)}</span>
      <input class="wslider foot-fixed" type="range" aria-label="Slab thickness">
      <input class="width num" type="number" min="0" step="any" aria-label="Slab thickness">`;
    shell.section.append(foot);
    const q = shell.q;
    const p = {
      fixed, x, y, step, lo, hi, version: 0, key, section: shell.section, zoom: null, drag: null, pinch: null, cutDrag: null, pointers: new Map(), more: [newLayer(), newLayer(), newLayer()],
      zoomReset: q('.zoom-reset'), slider: q('.slider'), center: q('.center'), wslider: q('.wslider'), width: q('.width'),
      caption: q('.view-caption'), pos: q('.view-pos'), angle: q('.view-angle'), canvas: q('canvas'), hover: q('.overlay-chip'),
    };
    p.canvas.setAttribute('aria-label', `${X.label}–${Y.label} intensity slice`);
    Object.assign(p.slider, { min: roundTo(lo, step), max: roundTo(hi, step), step, value: center });
    Object.assign(p.wslider, { min: step, max: roundTo(maxWidth, step), step, value: 3 * step });
    p.center.value = center;
    p.width.value = roundTo(3 * step, step);
    p.pos.textContent = `${F.label} = ${fmt(center)}`;
    p.slider.oninput = () => { p.center.value = p.slider.value; request(p); };
    p.center.oninput = () => { p.slider.value = p.center.value; paint(p.slider); request(p); };
    p.wslider.oninput = () => { p.width.value = p.wslider.value; request(p); };
    p.width.oninput = () => { p.wslider.value = p.width.value; paint(p.wslider); request(p); };
    // One pointer clicks or drags (by the click mode); two fingers pinch and pan.
    p.canvas.onpointerdown = (ev) => {
      if (ev.pointerType === 'touch') p.pointers.set(ev.pointerId, localPoint(p, ev));
      if (p.pointers.size === 2) startPinch(p, ev);
      else if (p.pointers.size < 2) (clickMode === 'cut' ? startCutDrag : startBox)(p, ev);
    };
    p.canvas.onpointermove = (ev) => {
      if (p.pointers.has(ev.pointerId)) p.pointers.set(ev.pointerId, localPoint(p, ev));
      if (p.pinch) movePinch(p);
      else {
        hover(p, ev);
        // A mouse released outside the page leaves no pointerup here.
        if (p.cutDrag && ev.pointerType === 'mouse' && !(ev.buttons & 1)) p.cutDrag = null;
        if (p.cutDrag) moveCutDrag(p, ev);
        else moveBox(p, ev);
      }
    };
    p.canvas.onpointerup = (ev) => {
      p.pointers.delete(ev.pointerId);
      // The finger left after a pinch does nothing until it lifts.
      if (p.pinch) { if (p.pointers.size < 2) p.pinch = null; }
      else if (p.cutDrag) p.cutDrag = null;
      else endBox(p, ev);
    };
    p.canvas.onpointercancel = (ev) => {
      p.pointers.delete(ev.pointerId);
      p.drag = p.pinch = p.cutDrag = null;
      p.canvas.classList.remove('panning');
      redraw();
    };
    p.canvas.onpointerleave = () => { p.hover.hidden = true; };
    p.canvas.ondblclick = () => { clearTimeout(p.clickTimer); setZoom(p, null); };
    p.zoomReset.onclick = () => setZoom(p, null);
    q('.save').onclick = () => savePNG(p);
    resized.observe(p.canvas);
    panels.push(p);
  }

  const shell = viewShell('3d', slotSwitch('3d'), 'Isosurface');
  shell.actions.innerHTML = `
    <button type="button" class="icon-btn reset3d" title="Reset camera">${ICONS.reset}</button>
    <button type="button" class="icon-btn" title="3-D options" aria-haspopup="dialog" aria-expanded="false" data-pop="pop-3d">${ICONS.gear}</button>
    <button type="button" class="icon-btn save" title="Save PNG">${ICONS.save}</button>`;
  shell.actions.append(shell.focusBtn, shell.maxBtn);
  shell.q('.view-body').innerHTML = '<canvas class="scene" role="img" aria-label="3-D isosurface with the three slices"></canvas>';
  const foot = document.createElement('footer');
  foot.className = 'view-foot';
  foot.innerHTML = `
    <span class="label">Level</span>
    <input class="iso-slider foot-grow" type="range" min="0" max="1000" aria-label="Isosurface level (log scale)">
    <input class="iso-level num" type="number" step="any" style="width: 7.6rem" aria-label="Isosurface level">
    <span class="foot-sep"></span>
    <span class="label">Surface</span><input class="iso-opacity foot-fixed short" type="range" min="0.05" max="1" step="0.05" value="0.6" aria-label="Surface opacity">
    <span class="label">Slices</span><input class="slice-opacity foot-fixed short" type="range" min="0.05" max="1" step="0.05" value="1" aria-label="Slice opacity">`;
  shell.section.append(foot);
  iso = {
    levelInput: shell.q('.iso-level'), slider: shell.q('.iso-slider'), opacity: shell.q('.iso-opacity'),
    sliceOpacity: shell.q('.slice-opacity'), caption: shell.q('.view-caption'), canvas: shell.q('canvas'),
    reset: shell.q('.reset3d'), save: shell.q('.save'), userLevel: false, range: null,
  };
  iso.caption.textContent = 'loading 3-D view…';
  setSegmented($('iso-grid'), '100');
  $('iso-slices').checked = true;
  setupPowder();
  setupCut();

  setLayout(layout, primary);
  paintAll();
}

/** The I(Q) view: a line plot of the powder average, with its shells set in the footer. */
function setupPowder() {
  const shell = viewShell('iq', slotSwitch('iq'), 'Powder average');
  shell.actions.innerHTML = `
    <button type="button" class="btn btn-ghost btn-xs zoom-reset" hidden title="Back to the full range (or double-click the plot)">Reset zoom</button>
    <div class="segmented sm pq-scale" role="group" aria-label="Intensity scale"><button type="button" data-value="linear" title="Linear intensity scale">Lin</button><button type="button" data-value="log" title="Logarithmic intensity scale">Log</button></div>
    <button type="button" class="icon-btn" title="I(Q) options" aria-haspopup="dialog" aria-expanded="false" data-pop="pop-iq">${ICONS.gear}</button>
    <button type="button" class="icon-btn data" title="Download I(Q) as text: Q, I, σ, coverage and voxels per shell">${ICONS.download}</button>
    <button type="button" class="icon-btn save" title="Save PNG">${ICONS.save}</button>`;
  shell.actions.append(shell.focusBtn, shell.maxBtn);
  shell.q('.view-body').innerHTML = `<canvas class="plot" role="img" aria-label="I(Q), the powder average of the volume"></canvas><span class="overlay-chip" hidden></span>
    <button type="button" class="btn btn-primary view-cta" title="Reduce the masked, symmetrized volume to I(Q)">Compute I(Q)</button>`;
  const foot = document.createElement('footer');
  foot.className = 'view-foot';
  const binsHelp = 'Q bins as Mantid Rebin parameters: a width in Å⁻¹ (0.05), a negative step for logarithmic bins (−0.01: ΔQ/Q = 1%), '
    + 'or ranges with their own steps (0.5, 0.02, 3, 0.05, 10). Empty for the shortest bin step in |Q|.';
  foot.innerHTML = `
    <span class="label" title="${binsHelp}">ΔQ</span>
    <input class="pq-slider foot-grow" type="range" min="0" max="1000" aria-label="Q bin width">
    <input class="pq-bins mono" type="text" spellcheck="false" autocomplete="off" aria-label="Q bins" title="${binsHelp}">
    <span class="foot-sep"></span>
    <span class="label">Q</span>
    <input class="pq-qmin num" type="number" min="0" step="any" aria-label="Q min" title="Smallest |Q| in Å⁻¹; empty for 0 (the shortest bin step for logarithmic bins)">
    <span class="unit">–</span>
    <input class="pq-qmax num" type="number" min="0" step="any" placeholder="all" aria-label="Q max" title="Largest |Q| in Å⁻¹; empty for all the data">
    <span class="unit">Å⁻¹</span>`;
  shell.section.append(foot);
  const q = shell.q;
  powder = {
    canvas: q('canvas'), hover: q('.overlay-chip'), caption: q('.view-caption'), zoomReset: q('.zoom-reset'),
    slider: q('.pq-slider'), bins: q('.pq-bins'), qmin: q('.pq-qmin'), qmax: q('.pq-qmax'), a: newCurveLayer(), more: [newCurveLayer(), newCurveLayer(), newCurveLayer()],
    compute: q('.view-cta'), requested: false, stale: true, download: false, plan: null, zoom: null, drag: null, view: null, at: null, step: null,
    // As a plot (see drawCurves()):
    xMin: 0, title: 'I(Q)',
    log: () => powderScale === 'log', band: () => $('iq-band').checked, cover: () => $('iq-coverage').checked,
    end: powderEnd, xLabel: () => 'Q (Å⁻¹)', yLabel: (log) => (log ? 'I(Q), log scale' : 'I(Q)'),
    empty: () => {
      // Before it is asked for, the view holds only its Compute button.
      if (!powder.compute.hidden) return null;
      const busy = [powder.a, ...powder.more].some((l) => l.busy), failed = layerOf(powder, shownIndex()).error;
      return { text: busy ? 'Computing I(Q)…' : failed || 'I(Q) appears here once computed.', error: !busy && !!failed };
    },
    extras: (data) => [binsLabel(powder.plan.bins), `voxels split ${data.split}³`, data.order > 1 ? data.symmetry : 'no symmetry',
      data.masked && mask ? `mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : 'no mask'],
    // The shell: Q, d = 2π/Q, and I ± σ with the coverage for each dataset.
    readout: (b, layers) => {
      const value = (d) => (Number.isFinite(d.intensity[b])
        ? `${fmtValue(d.intensity[b])}${Number.isFinite(d.sigma[b]) ? ` ± ${fmtValue(d.sigma[b])}` : ''} (${pct(d.coverage[b])} covered)`
        : 'no data');
      const Q = shellMid(layers[0].data, b);
      return `Q ${fmt(Q)} Å⁻¹  d ${fmt(2 * Math.PI / Q, 3)} Å  →  ${layers.map((l) => `${layers.length > 1 ? `${l.letter} ` : ''}${value(l.data)}`).join(' · ')}`;
    },
  };
  powder.compute.onclick = () => {
    powder.requested = true;
    showPowder();
    updatePowder();
    drawPlot(powder);
  };
  try {
    powder.step = qExtent(meta.dims, meta.lattice).step;
    powder.bins.placeholder = String(Number(powder.step.toPrecision(1)));
  } catch (err) {
    powder.a.error = err.message;
  }
  syncBinsSlider();
  setSegmented(q('.pq-scale'), powderScale);
  segmented(q('.pq-scale'), (value) => {
    powderScale = value;
    // The intensity axis changes, so a zoom keeps its Q range only.
    if (powder.zoom) setCurveZoom(powder, powder.zoom.x ? { x: powder.zoom.x, y: null } : null);
    persist();
    drawPlot(powder);
  });
  // The slider recomputes as it moves (one request in flight, the latest one next); typed values apply on Enter or leaving the field.
  powder.slider.oninput = slideBins;
  powder.bins.onchange = () => { syncBinsSlider(); requestPowder(); };
  powder.qmin.onchange = powder.qmax.onchange = () => requestPowder();
  for (const input of [powder.bins, powder.qmin, powder.qmax]) input.onkeydown = (e) => { if (e.key === 'Enter') input.blur(); };
  wirePlot(powder);
  q('.data').onclick = downloadPowder;
  q('.save').onclick = () => savePlotPNG(powder, `${shownStem()}_IQ.png`);
  showPowder();
}

// ---- Layout ------------------------------------------------------------------------

// Below this width the views stack, before the panel, and the layouts do not apply.
const compactLayout = matchMedia('(max-width: 1000px)');
compactLayout.addEventListener('change', () => { if (views.hk) setLayout(layout, primary); });

function setLayout(mode, key = primary) {
  if (!views[key]) key = 'hk';
  if (mode !== 'single') lastMulti = mode;
  if (SLOT_VIEWS.includes(key)) slot = key;
  layout = mode;
  primary = key;
  const ws = $('workspace');
  ws.classList.remove('quad', 'focus', 'single');
  ws.classList.add(mode);
  setSegmented($('layout'), mode);
  for (const [k, v] of Object.entries(views)) {
    const isPrimary = k === key;
    v.section.classList.toggle('primary', isPrimary);
    // Small views in the focus layout are thumbnails (styled by .thumb), reachable with Tab and opened with Enter.
    const thumb = mode === 'focus' && !isPrimary && !compactLayout.matches;
    v.section.classList.toggle('thumb', thumb);
    if (thumb) v.section.tabIndex = 0;
    else v.section.removeAttribute('tabindex');
    const back = mode === 'focus' && isPrimary;
    v.focusBtn.innerHTML = back ? ICONS.quad : ICONS.focus;
    v.focusBtn.title = back ? 'Back to four views' : 'Show this view large';
    const restore = mode === 'single' && isPrimary;
    v.maxBtn.innerHTML = restore ? ICONS.restore : ICONS.max;
    v.maxBtn.title = restore ? 'Restore the layout (Esc)' : 'Maximize this view (double-click the header)';
    v.section.classList.toggle('off-slot', SLOT_VIEWS.includes(k) && k !== slot);
  }
  try { localStorage.setItem('nxv-layout', JSON.stringify({ mode: lastMulti, key, slot })); } catch { /* storage unavailable */ }
  updatePowder();
  updateCut();
}

/** Show the 3-D view, I(Q) or the cut (`key`) in the fourth place. */
function setSlot(key) {
  slot = key;
  setLayout(layout, SLOT_VIEWS.includes(primary) ? key : primary);
}

// Redraw a plot whenever its canvas changes size (layout switches, window resizes),
// and compute I(Q) once its view is shown.
let resizeQueued = false;
const resized = new ResizeObserver(() => {
  if (resizeQueued) return;
  resizeQueued = true;
  requestAnimationFrame(() => { resizeQueued = false; redraw(); updatePowder(); updateCut(); });
});

// ---- Box zoom and clicks -------------------------------------------------------------

/** The pointer in the plot's drawing coordinates (CSS px over the UI scale). */
function localPoint(p, e) {
  const box = p.canvas.getBoundingClientRect(), s = uiScale();
  return [(e.clientX - box.left) / s, (e.clientY - box.top) / s];
}

function startBox(p, e) {
  if (e.button !== 0 || !p.inverse) return;
  const [x, y] = localPoint(p, e);
  // Keep the mapping at the start of the drag: panning changes p.inverse as it goes.
  p.drag = { x0: x, y0: y, x1: x, y1: y, moved: false, view: p.view, inverse: p.inverse };
  try { p.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic or already released pointer */ }
  if (clickMode === 'move') p.canvas.classList.add('panning');
}

function setClickMode(mode) {
  clickMode = mode;
  setSegmented($('click-mode'), mode);
  $('workspace').dataset.click = mode;
  // The cut shows on its slice in Cut mode, and its empty view says what to do.
  if (!cut) return;
  showCut();
  drawPlot(cut);
  const p = cutPanel();
  if (p) drawPanel(p);
}

// Navigate mode: dragging keeps moving the other slices. Zoom mode: dragging
// draws the zoom box. Move mode: dragging slides the visible region.
function moveBox(p, e) {
  if (!p.drag) return;
  [p.drag.x1, p.drag.y1] = localPoint(p, e);
  p.drag.moved ||= Math.hypot(p.drag.x1 - p.drag.x0, p.drag.y1 - p.drag.y0) > (clickMode === 'move' ? 1 : 5);
  if (!p.drag.moved) return;
  if (clickMode === 'zoom') drawPanel(p);
  else if (clickMode === 'move') pan(p, p.drag);
  else navigate(p, e);
}

/**
 * Shift the view window by the dragged distance, measured with the mapping
 * from the start of the drag. The window center stays inside the data.
 */
function pan(p, d) {
  const [a0, b0] = d.inverse(d.x0, d.y0), [a1, b1] = d.inverse(d.x1, d.y1);
  const [u0, u1, v0, v1] = d.view, X = meta.dims[p.x], Y = meta.dims[p.y];
  const shift = (lo, hi, delta, edges) => {
    const c = clamp((lo + hi) / 2 - delta, edges[0], edges[edges.length - 1]);
    return [c - (hi - lo) / 2, c + (hi - lo) / 2];
  };
  setZoom(p, { u: shift(u0, u1, a1 - a0, X.edges), v: shift(v0, v1, b1 - b0, Y.edges) });
}

function endBox(p, e) {
  const d = p.drag;
  p.drag = null;
  p.canvas.classList.remove('panning');
  if (!d) return;
  if (clickMode === 'move') { if (d.moved) pan(p, d); return; }
  if (d.moved) {
    if (clickMode === 'navigate') { navigate(p, e); return; }
    const box = boxWindow(p, d);
    if (box) setZoom(p, box);
    else drawPanel(p);
    return;
  }
  // A click: wait briefly so a double-click (full view) does not also act on it.
  const { clientX, clientY } = e;
  clearTimeout(p.clickTimer);
  p.clickTimer = setTimeout(() => (clickMode === 'zoom' ? zoomAt : navigate)(p, { clientX, clientY }), 250);
}

/** A second finger: drop the click or drag of the first, and pinch from here. */
function startPinch(p, e) {
  if (!p.inverse) return;
  clearTimeout(p.clickTimer);
  p.drag = p.cutDrag = null;
  p.canvas.classList.remove('panning');
  try { p.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic or already released pointer */ }
  const [a, b] = [...p.pointers.values()];
  p.pinch = { a, b, view: p.view, inverse: p.inverse };
  drawPanel(p);
}

/**
 * Scale the window by the change in finger spread and keep the data point that
 * was under their midpoint under it, with the mapping from the start of the pinch:
 * the window center maps to the plot center, so a window k times larger, centered
 * at m - k (inverse(mid) - center), puts m under mid.
 */
function movePinch(p) {
  const d = p.pinch, [a, b] = [...p.pointers.values()];
  const spread = (s, t) => Math.hypot(t[0] - s[0], t[1] - s[1]);
  const middle = (s, t) => [(s[0] + t[0]) / 2, (s[1] + t[1]) / 2];
  if (!b || spread(a, b) < 1) return;
  const [u0, u1, v0, v1] = d.view, X = meta.dims[p.x], Y = meta.dims[p.y];
  const ex = [X.edges[0], X.edges[X.edges.length - 1]], ey = [Y.edges[0], Y.edges[Y.edges.length - 1]];
  // At least two bins across, and at most the whole grid, on both axes.
  const least = Math.max(2 * (ex[1] - ex[0]) / (X.edges.length - 1) / (u1 - u0), 2 * (ey[1] - ey[0]) / (Y.edges.length - 1) / (v1 - v0));
  const most = Math.min((ex[1] - ex[0]) / (u1 - u0), (ey[1] - ey[0]) / (v1 - v0));
  const k = clamp(spread(d.a, d.b) / spread(a, b), Math.min(least, most), most);
  const [mu, mv] = d.inverse(...middle(d.a, d.b)), [nu, nv] = d.inverse(...middle(a, b));
  const cu = mu - k * (nu - (u0 + u1) / 2), cv = mv - k * (nv - (v0 + v1) / 2);
  const hu = k * (u1 - u0) / 2, hv = k * (v1 - v0) / 2;
  const fit = (c, h, [lo, hi]) => { const m = clamp(c, lo, hi); return [m - h, m + h]; };
  setZoom(p, { u: fit(cu, hu, ex), v: fit(cv, hv, ey) });
}

/** Zoom in by `factor` around the clicked point, keeping the window inside the data. */
function zoomAt(p, e, factor = 2) {
  const pt = pointAt(p, e);
  if (!pt) return;
  const [u0, u1, v0, v1] = p.view, X = meta.dims[p.x], Y = meta.dims[p.y];
  const fit = (c, half, edges) => {
    const lo = edges[0], hi = edges[edges.length - 1];
    if (2 * half >= hi - lo) return [lo, hi];
    const start = clamp(c - half, lo, hi - 2 * half);
    return [start, start + 2 * half];
  };
  const u = fit(pt[0], (u1 - u0) / (2 * factor), X.edges), v = fit(pt[1], (v1 - v0) / (2 * factor), Y.edges);
  // Stop at two bins across.
  const wx = 2 * (X.edges[X.edges.length - 1] - X.edges[0]) / (X.edges.length - 1);
  const wy = 2 * (Y.edges[Y.edges.length - 1] - Y.edges[0]) / (Y.edges.length - 1);
  if (u[1] - u[0] >= wx && v[1] - v[0] >= wy) setZoom(p, { u, v });
}

/** The (u, v) window enclosing a dragged screen rectangle, clipped to the data. */
function boxWindow(p, d) {
  const corners = [[d.x0, d.y0], [d.x1, d.y0], [d.x1, d.y1], [d.x0, d.y1]].map(([x, y]) => p.inverse(x, y));
  const X = meta.dims[p.x], Y = meta.dims[p.y];
  const ex = [X.edges[0], X.edges[X.edges.length - 1]], ey = [Y.edges[0], Y.edges[Y.edges.length - 1]];
  const u = [Math.max(ex[0], Math.min(...corners.map((c) => c[0]))), Math.min(ex[1], Math.max(...corners.map((c) => c[0])))];
  const v = [Math.max(ey[0], Math.min(...corners.map((c) => c[1]))), Math.min(ey[1], Math.max(...corners.map((c) => c[1])))];
  // At least two bins across each axis.
  const wx = 2 * (ex[1] - ex[0]) / (X.edges.length - 1), wy = 2 * (ey[1] - ey[0]) / (Y.edges.length - 1);
  if (!(u[1] - u[0] >= wx && v[1] - v[0] >= wy)) return null;
  // With a common scale, widen the shorter side so both axes span the same
  // length (an equal-sided window), centered on the dragged box.
  const g = geometry(p);
  if (g.equal) {
    const lu = (u[1] - u[0]) * g.lx, lv = (v[1] - v[0]) * g.ly;
    const widen = (r, len, l) => { const c = (r[0] + r[1]) / 2, h = len / l / 2; return [c - h, c + h]; };
    if (lu < lv) return { u: widen(u, lv, g.lx), v };
    return { u, v: widen(v, lu, g.ly) };
  }
  return { u, v };
}

function setZoom(p, box) {
  p.zoom = box;
  p.zoomReset.hidden = !box;
  drawPanel(p);
}

// ---- Color scale -------------------------------------------------------------

function readSettings() {
  const limit = $('limit').value.trim();
  const s = {
    cmap: $('cmap').value, scale: $('scale').dataset.value, min: Number($('vmin').value), max: Number($('vmax').value),
    soft: Number($('soft').value), limit: limit === '' ? Infinity : Number(limit),
    angles: $('angles').checked ? 'nominal' : 'measured', guides: $('guides').checked, grid: $('grid').checked,
  };
  if (![s.min, s.max].every(Number.isFinite) || s.max <= s.min) throw new Error('Use finite color limits with vmax > vmin.');
  if (s.scale === 'asinh' && !(s.soft > 0)) throw new Error('Asinh softening must be positive.');
  if (s.scale === 'log' && s.min <= 0) throw new Error('Log scale needs vmin > 0.');
  if (!(s.limit > 0)) throw new Error('View range must be positive; leave it empty for the full range.');
  return s;
}

function scaler(s) {
  const f = s.scale === 'asinh' ? (v) => Math.asinh(v / s.soft) : s.scale === 'log' ? Math.log10 : (v) => v;
  const lo = f(s.min), span = f(s.max) - lo;
  return (v) => {
    if (s.scale === 'log' && v <= 0) return 0;
    const t = (f(v) - lo) / span;
    return t < 0 ? 0 : t > 1 ? 1 : t;
  };
}

function autoRange() {
  // Signed data (a ΔPDF) uses the magnitudes of all values, on a range symmetric about 0.
  const positive = [], signed = !!meta.signed;
  for (const p of panels) {
    for (const layer of [p, ...readyOthers().map((o) => p.more[o.slot])]) {
      for (const v of layer.data?.values ?? []) if (v > 0 || (signed && v < 0)) positive.push(Math.abs(v));
    }
  }
  positive.sort((a, b) => a - b);
  const quantile = (q) => positive[Math.min(positive.length - 1, Math.floor(q * positive.length))];
  // Bragg peaks dominate the top percentiles, so a 97th-percentile ceiling
  // with the median as asinh softening keeps diffuse intensity visible.
  const vmax = positive.length ? sig(quantile(0.97)) : 1;
  $('vmin').value = $('scale').dataset.value === 'log' ? sig(vmax / 1000) : signed ? -vmax : 0;
  $('vmax').value = vmax;
  $('soft').value = positive.length ? sig(quantile(0.5), 1) : sig(vmax / 20);
  rangeIsAuto = true;
  redraw();
}

function colorTicks(s) {
  const lo = s.min, hi = s.max;
  if (s.scale === 'linear') return niceTicks(lo, hi, 4);
  const out = [];
  if (s.scale === 'asinh' && lo <= 0 && hi >= 0) out.push(0);
  const top = Math.ceil(Math.log10(Math.max(Math.abs(lo), Math.abs(hi))));
  const bottom = s.scale === 'log' ? Math.floor(Math.log10(lo)) : Math.floor(Math.log10(s.soft));
  for (let k = bottom; k <= top; k++) for (const v of [10 ** k, -(10 ** k)]) if (v >= lo && v <= hi) out.push(v);
  return out.sort((a, b) => a - b);
}

function paintColorbar() {
  const lut = LUTS[$('cmap').value];
  for (const id of ['legend-bar']) {
    const c = $(id).getContext('2d'), img = c.createImageData(256, 1);
    for (let k = 0; k < 256; k++) img.data.set([lut[3 * k], lut[3 * k + 1], lut[3 * k + 2], 255], 4 * k);
    c.putImageData(img, 0, 0);
  }
}

// ---- Drawing -------------------------------------------------------------------

function geometry(p) {
  const cell = meta.lattice && (settings.angles === 'nominal' ? nominalCell(meta.lattice) : meta.lattice);
  return planeGeometry(meta.dims, cell, p.x, p.y);
}

function viewRange(dim, limit) {
  const e = dim.edges;
  const lo = Math.max(e[0], -limit), hi = Math.min(e[e.length - 1], limit);
  return lo < hi ? [lo, hi] : [e[0], e[e.length - 1]];
}

/** The colored image of one dataset's slice (a panel, or another dataset's layer in `p.more`), cached. */
function layerImage(p, s) {
  const key = `${p.version}|${s.cmap}|${s.scale}|${s.min}|${s.max}|${s.soft}`;
  if (p.imageKey === key) return p.image;
  const { values, rows, cols } = p.data;
  const image = p.image ?? document.createElement('canvas');
  image.width = cols;
  image.height = rows;
  const ctx = image.getContext('2d');
  const pixels = ctx.createImageData(cols, rows), d = pixels.data, lut = LUTS[s.cmap], norm = scaler(s);
  for (let i = 0, o = 0; i < values.length; i++, o += 4) {
    const v = values[i];
    if (v === v) {
      const k = Math.round(norm(v) * 255) * 3;
      d[o] = lut[k]; d[o + 1] = lut[k + 1]; d[o + 2] = lut[k + 2];
      d[o + 3] = 255;
    } else {
      // Transparent: the 2-D plot shows its grey fill, the 3-D view a hole.
      d[o + 3] = 0;
    }
  }
  ctx.putImageData(pixels, 0, 0);
  p.image = image;
  p.imageKey = key;
  return image;
}

/**
 * Draw a panel into `canvas` (w x h CSS px). On screen the shared colorbar
 * lives in the toolbar; exports add a title and their own colorbar, without guides.
 */
function draw(p, canvas, w, h, dpr, exporting = false) {
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const c = canvas.getContext('2d');
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  c.fillStyle = '#ffffff';
  c.fillRect(0, 0, w, h);
  if (!p.data || !settings || w < 40 || h < 40) return;

  const s = settings, X = meta.dims[p.x], Y = meta.dims[p.y], F = meta.dims[p.fixed];
  const g = geometry(p), sin = Math.sqrt(Math.max(0, 1 - g.cos * g.cos));
  const [u0, u1] = p.zoom?.u ?? viewRange(X, s.limit), [v0, v1] = p.zoom?.v ?? viewRange(Y, s.limit);
  const wx = (u, v) => g.lx * u + g.ly * g.cos * v, wy = (v) => g.ly * sin * v;
  const box = [[u0, v0], [u1, v0], [u1, v1], [u0, v1]];
  const xs = box.map(([u, v]) => wx(u, v)), ys = box.map(([, v]) => wy(v));
  const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys);
  const pad = { l: 58, r: exporting ? 76 : 20, t: exporting ? 40 : 14, b: 46 };
  const aw = Math.max(10, w - pad.l - pad.r), ah = Math.max(10, h - pad.t - pad.b);
  let sx = aw / (xmax - xmin), sy = ah / (ymax - ymin);
  if (g.equal) sx = sy = Math.min(sx, sy);
  const ox = pad.l + (aw - sx * (xmax - xmin)) / 2 - sx * xmin;
  const oy = pad.t + (ah - sy * (ymax - ymin)) / 2 + sy * ymax;
  const project = (u, v) => [ox + sx * wx(u, v), oy - sy * wy(v)];
  if (!exporting) {
    p.view = [u0, u1, v0, v1];
    p.project = project;
    p.inverse = (px, py) => {
      const v = (oy - py) / (sy * g.ly * sin);
      return [((px - ox) / sx - g.ly * g.cos * v) / g.lx, v];
    };
  }

  // Intensity images, clipped to the view parallelogram, or to each dataset's
  // part of it when comparing (see splitLayout()).
  const corners = box.map(([u, v]) => project(u, v));
  const path = (...polygons) => {
    c.beginPath();
    for (const pts of polygons) {
      pts.forEach(([a, b], i) => (i ? c.lineTo(a, b) : c.moveTo(a, b)));
      c.closePath();
    }
  };
  path(corners);
  c.fillStyle = MISSING;
  c.fill();
  const shown = compareShown(), split = splitLayout(u0, u1, v0, v1), onCanvas = (pts) => pts.map(([u, v]) => project(u, v));
  if (split.hatch.length) {
    // The other datasets' parts are hatched along the diagonal, so they stand out where they have no data.
    c.save();
    path(...split.hatch.map(onCanvas));
    c.clip();
    const [ax, ay] = corners[3], len = Math.hypot(corners[1][0] - ax, corners[1][1] - ay);
    const d = [(corners[1][0] - ax) / len, (corners[1][1] - ay) / len], n = [-d[1], d[0]];
    const offsets = corners.map(([x, y]) => (x - ax) * n[0] + (y - ay) * n[1]);
    c.strokeStyle = 'rgba(71, 84, 103, 0.3)';
    c.lineWidth = 1;
    c.beginPath();
    for (let t = Math.min(...offsets); t <= Math.max(...offsets); t += 6) {
      const [x, y] = [ax + t * n[0], ay + t * n[1]];
      c.moveTo(x - 2 * len * d[0], y - 2 * len * d[1]);
      c.lineTo(x + 2 * len * d[0], y + 2 * len * d[1]);
    }
    c.stroke();
    c.restore();
  }
  for (const [k, part] of split.parts) {
    const layer = layerOf(p, k), dims = k ? slots[k - 1].meta.dims : meta.dims;
    if (!layer.data) continue;
    const { rows, cols } = layer.data, ex = dims[p.x].edges, ey = dims[p.y].edges;
    const dx = (ex[ex.length - 1] - ex[0]) / cols, dy = (ey[ey.length - 1] - ey[0]) / rows;
    c.save();
    path(onCanvas(part));
    c.clip();
    const [x0, y0] = project(ex[0], ey[0]);
    c.transform(sx * g.lx * dx, 0, sx * g.ly * g.cos * dy, -sy * g.ly * sin * dy, x0, y0);
    c.imageSmoothingEnabled = false;
    c.drawImage(layerImage(layer, s), 0, 0);
    c.restore();
  }

  // Integer grid: lines at whole-number values of both axes (integer H, K, L on
  // r.l.u. axes), along the true axis directions. An axis whose lines would be
  // closer than 6 px apart on screen gets none.
  if (s.grid) {
    const gapU = sx * g.lx * sy * sin / Math.hypot(sx * g.cos, sy * sin), gapV = sy * g.ly * sin;
    const lines = [];
    if (gapU >= 6) for (let n = Math.ceil(u0); n <= u1; n++) lines.push([project(n, v0), project(n, v1)]);
    if (gapV >= 6) for (let n = Math.ceil(v0); n <= v1; n++) lines.push([project(u0, n), project(u1, n)]);
    c.save();
    path(corners);
    c.clip();
    // A light line on a faint dark one, so the grid shows on dark and bright colors alike.
    for (const [color, width] of [['rgba(18, 24, 33, 0.28)', 2], ['rgba(255, 255, 255, 0.55)', 0.75]]) {
      c.strokeStyle = color;
      c.lineWidth = width;
      c.beginPath();
      for (const [[a, b], [e, f]] of lines) { c.moveTo(a, b); c.lineTo(e, f); }
      c.stroke();
    }
    c.restore();
  }

  // Dashed lines where the other two slices cut this plane.
  if (s.guides && !exporting) {
    const at = (dim) => panels.find((q) => q.fixed === dim)?.data?.center;
    const gu = at(p.x), gv = at(p.y), lines = [];
    if (gu >= u0 && gu <= u1) lines.push([project(gu, v0), project(gu, v1)]);
    if (gv >= v0 && gv <= v1) lines.push([project(u0, gv), project(u1, gv)]);
    c.save();
    c.lineWidth = 1;
    c.setLineDash([4, 4]);
    for (const [[a, b], [e, f]] of lines) {
      for (const [color, offset] of [['rgba(255,255,255,.9)', 0], ['rgba(18,24,33,.55)', 4]]) {
        c.strokeStyle = color;
        c.lineDashOffset = offset;
        c.beginPath(); c.moveTo(a, b); c.lineTo(e, f); c.stroke();
      }
    }
    c.restore();
  }

  // The cuts between the datasets' parts (a white gap with dark edges) and dataset tags when comparing.
  if (split.lines.length) {
    c.save();
    c.lineCap = 'butt';
    c.shadowColor = 'rgba(18, 24, 33, 0.25)';
    c.shadowBlur = 3;
    for (const [color, width] of [['rgba(18, 24, 33, 0.75)', 4.5], ['#ffffff', 2.5]]) {
      c.strokeStyle = color;
      c.lineWidth = width;
      c.beginPath();
      for (const [from, to] of split.lines) { c.moveTo(...project(...from)); c.lineTo(...project(...to)); }
      c.stroke();
      c.shadowColor = 'transparent';
    }
    c.restore();
  }
  if (shown) {
    // Each tag sits in its part's corner of the view: [corner, along, side, left, down].
    const room = 0.45 * Math.hypot(corners[1][0] - corners[0][0], corners[1][1] - corners[0][1]);
    const at = { 0: [0, 1, 3, false, false], 1: [1, 0, 2, true, false], 2: [2, 3, 1, true, true], 3: [3, 2, 0, false, true] };
    for (const [k, , corner] of split.parts) {
      const [i, along, side, left, down] = at[corner];
      datasetTag(c, LETTERS[k], k ? slots[k - 1].name : sourceName, corners[i], corners[along], corners[side], room, left, down);
    }
  }
  if (!exporting) drawCutOverlay(c, p, corners);

  // Axes, ticks and labels.
  c.strokeStyle = AXIS;
  c.lineWidth = 1;
  c.beginPath();
  c.moveTo(...corners[3]); c.lineTo(...corners[0]); c.lineTo(...corners[1]);
  c.stroke();
  c.fillStyle = INK2;
  c.font = `10.5px ${MONO}`;
  c.textAlign = 'center';
  c.textBaseline = 'top';
  // Axes that share a scale share a tick step, so equal lengths look equal.
  const step = g.equal ? Math.max(niceStep(Math.min(u1 - u0, v1 - v0) / 5), niceStep(Math.max(u1 - u0, v1 - v0) / 10)) : null;
  for (const t of niceTicks(u0, u1, 5, step)) {
    const [a, b] = project(t, v0);
    c.beginPath(); c.moveTo(a, b); c.lineTo(a, b + 4); c.stroke();
    c.fillText(fmt(t), a, b + 7);
  }
  c.textAlign = 'right';
  c.textBaseline = 'middle';
  for (const t of niceTicks(v0, v1, 5, step)) {
    const [a, b] = project(u0, t);
    c.beginPath(); c.moveTo(a, b); c.lineTo(a - 4, b); c.stroke();
    c.fillText(fmt(t), a - 7, b);
  }
  c.fillStyle = INK;
  c.font = `600 11.5px ${SANS}`;
  c.textAlign = 'center';
  const [bx, by] = project((u0 + u1) / 2, v0);
  c.fillText(withUnits(X), bx, by + 31);
  const [ax0, ay0] = project(u0, v0), [ax1, ay1] = project(u0, v1);
  const len = Math.hypot(ax1 - ax0, ay1 - ay0), nx = (ay1 - ay0) / len, ny = -(ax1 - ax0) / len;
  c.save();
  c.translate((ax0 + ax1) / 2 + 44 * nx, (ay0 + ay1) / 2 + 44 * ny);
  c.rotate(Math.atan2(ay1 - ay0, ax1 - ax0));
  c.fillText(withUnits(Y), 0, 0);
  c.restore();

  if (!exporting) return;
  c.textAlign = 'left';
  c.textBaseline = 'alphabetic';
  c.font = `650 14px ${SANS}`;
  const d = layerOf(p, shownIndex()).data ?? p.data;
  const title = `${F.label} = ${fmt(d.center)}`;
  c.fillText(title, 10, 24);
  const titleWidth = c.measureText(title).width;
  c.font = `11px ${MONO}`;
  c.fillStyle = INK2;
  const extras = `${symmetryNote(d)}${mask ? ` · mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : ''}${d.removed ? ' · removed only' : ''}`;
  c.fillText(`slab ${fmt(d.slab[0])} to ${fmt(d.slab[1])}${extras}`, 20 + titleWidth, 24);

  // Colorbar (exports only).
  const lut = LUTS[s.cmap], cx = w - pad.r + 22, cw = 10, top = pad.t + 4, bottom = h - pad.b, span = bottom - top;
  for (let k = 0; k < 256; k++) {
    c.fillStyle = `rgb(${lut[3 * k]},${lut[3 * k + 1]},${lut[3 * k + 2]})`;
    c.fillRect(cx, bottom - (k + 1) * span / 256, cw, span / 256 + 0.6);
  }
  c.strokeStyle = AXIS;
  c.strokeRect(cx, top, cw, span);
  c.fillStyle = INK2;
  c.font = `10px ${MONO}`;
  c.textAlign = 'left';
  c.textBaseline = 'middle';
  const norm = scaler(s), yOf = (v) => bottom - norm(v) * span;
  const labels = [s.min, s.max];
  for (const v of colorTicks(s)) if (labels.every((k) => Math.abs(yOf(k) - yOf(v)) >= 14)) labels.push(v);
  for (const v of labels) {
    const y = yOf(v);
    c.beginPath(); c.moveTo(cx + cw, y); c.lineTo(cx + cw + 3, y); c.stroke();
    c.fillText(fmtValue(v), cx + cw + 6, y);
  }
}

/**
 * How the visible window [u0, u1] × [v0, v1] is shared when comparing, in plot
 * coordinates: `parts`, [k, polygon, corner] for each dataset k shown (0 for A)
 * with the view corner (0–3, from (u0, v0) counterclockwise) where its tag
 * goes; `lines` that divide the parts; and `hatch`, the polygons outside A's
 * part. Two datasets split the window along the diagonal from (u0, v1) to
 * (u1, v0), A below; three or four take its quadrants about the center, A lower
 * left, B lower right, C upper left and D upper right (a quadrant stays empty
 * until its dataset is loaded). One dataset alone fills it.
 */
function splitLayout(u0, u1, v0, v1) {
  const shown = compareShown(), box = [[u0, v0], [u1, v0], [u1, v1], [u0, v1]];
  if (shown !== 'split') {
    const k = shownIndex();
    return { parts: [[k, box, k ? 2 : 0]], lines: [], hatch: k ? [box] : [] };
  }
  const ready = [0, ...readyOthers().map((d) => d.slot + 1)];
  if (others().length < 2) {
    const above = [box[1], box[2], box[3]];
    return { parts: [[0, [box[0], box[1], box[3]], 0], [ready[1], above, 2]], lines: [[box[3], box[1]]], hatch: [above] };
  }
  const um = (u0 + u1) / 2, vm = (v0 + v1) / 2;
  const quad = (a, b, c, d) => [[a, c], [b, c], [b, d], [a, d]];
  const quadrants = [quad(u0, um, v0, vm), quad(um, u1, v0, vm), quad(u0, um, vm, v1), quad(um, u1, vm, v1)];
  return {
    parts: ready.map((k) => [k, quadrants[k], [0, 1, 3, 2][k]]),
    lines: [[[um, v0], [um, v1]], [[u0, vm], [u1, vm]]],
    hatch: quadrants.slice(1),
  };
}

/**
 * A dataset tag (letter and file name) inside the view corner `at`, offset
 * toward its neighbouring corners `along` (same row) and `side` (same column).
 * The tag grows right and up from the corner, or to the left / down when `left` / `down`.
 */
function datasetTag(c, letter, name, at, along, side, room, left = false, down = false) {
  const unit = ([x, y]) => {
    const n = Math.hypot(x - at[0], y - at[1]) || 1;
    return [(x - at[0]) / n, (y - at[1]) / n];
  };
  const [ax, ay] = unit(along), [bx, by] = unit(side);
  const x = at[0] + 12 * ax + 8 * bx, y = at[1] + 12 * ay + 8 * by;
  c.save();
  c.font = `600 12px ${SANS}`;
  const h = 26, badge = 20, maxText = room - badge - 16;
  // A long name keeps the part that differs from another dataset's name (A's, or B's for A).
  const fits = (t) => c.measureText(t).width <= maxText;
  const text = maxText > 24 ? shortenName(name, fits, name === sourceName ? others()[0]?.name : sourceName) : '';
  const tw = text && fits(text) ? c.measureText(text).width : 0;
  const w = 3 + badge + (tw ? 7 + tw + 9 : 3), x0 = left ? x - w : x, top = down ? y : y - h;
  c.fillStyle = 'rgba(255, 255, 255, 0.9)';
  c.strokeStyle = 'rgba(18, 24, 33, 0.18)';
  c.lineWidth = 1;
  c.beginPath(); c.roundRect(x0, top, w, h, 6); c.fill(); c.stroke();
  c.fillStyle = INK;
  c.beginPath(); c.roundRect(x0 + 3, top + 3, badge, badge, 5); c.fill();
  c.fillStyle = '#ffffff';
  c.textAlign = 'center';
  c.textBaseline = 'middle';
  c.font = `700 12.5px ${SANS}`;
  c.fillText(letter, x0 + 3 + badge / 2, top + h / 2 + 0.5);
  if (tw) {
    c.font = `600 12px ${SANS}`;
    c.fillStyle = INK2;
    c.textAlign = 'left';
    c.fillText(text, x0 + 3 + badge + 7, top + h / 2 + 0.5);
  }
  c.restore();
}

/**
 * The UI scale: 1 rem over 10 CSS px, above 1 on large screens at 100% (set in
 * index.html). Plots are drawn at this scale, so their text grows with the page's.
 */
const uiScale = () => parseFloat(getComputedStyle(document.documentElement).fontSize) / 10 || 1;

/** Redraw one panel, with the zoom box being dragged (as its u-v parallelogram). */
function drawPanel(p) {
  if (!settings) return;
  const s = uiScale();
  draw(p, p.canvas, p.canvas.clientWidth / s, p.canvas.clientHeight / s, s * (devicePixelRatio || 1));
  const d = p.drag;
  if (!d?.moved || !p.project || clickMode !== 'zoom') return;
  const corners = [[d.x0, d.y0], [d.x1, d.y0], [d.x1, d.y1], [d.x0, d.y1]].map(([x, y]) => p.inverse(x, y));
  const u = [Math.min(...corners.map((c) => c[0])), Math.max(...corners.map((c) => c[0]))];
  const v = [Math.min(...corners.map((c) => c[1])), Math.max(...corners.map((c) => c[1]))];
  const c = p.canvas.getContext('2d');
  c.save();
  c.beginPath();
  [[u[0], v[0]], [u[1], v[0]], [u[1], v[1]], [u[0], v[1]]].map(([a, b]) => p.project(a, b)).forEach(([x, y], i) => (i ? c.lineTo(x, y) : c.moveTo(x, y)));
  c.closePath();
  c.fillStyle = 'rgba(47, 116, 230, 0.12)';
  c.fill();
  c.setLineDash([5, 3]);
  c.strokeStyle = '#2f74e6';
  c.lineWidth = 1.5;
  c.stroke();
  c.restore();
}

function redraw() {
  if (!meta || !panels.length) return;
  drawPlot(powder);
  drawPlot(cut);
  try {
    settings = readSettings();
    error('');
  } catch (err) {
    error(err.message);
    return;
  }
  $('soft-wrap').hidden = settings.scale !== 'asinh';
  for (const p of panels) {
    const g = geometry(p), angle = Math.acos(clamp(g.cos, -1, 1)) * 180 / Math.PI;
    p.angle.hidden = !g.lattice || Math.abs(angle - 90) < 0.5;
    p.angle.textContent = `∠ ${angle.toFixed(1)}°`;
    p.angle.title = `Angle between the ${meta.dims[p.x].label} and ${meta.dims[p.y].label} axes (${isDirect(meta.dims[p.x]) ? 'direct' : 'reciprocal'} lattice, ${settings.angles} cell angles)`;
    drawPanel(p);
  }
  paintColorbar();
  updateStates();
  update3D();
}

// ---- Interaction -----------------------------------------------------------------

function pointAt(p, e) {
  if (!p.data || !p.inverse) return null;
  const [u, v] = p.inverse(...localPoint(p, e));
  const [u0, u1, v0, v1] = p.view;
  return u >= u0 && u <= u1 && v >= v0 && v <= v1 ? [u, v] : null;
}

function hover(p, e) {
  const pt = pointAt(p, e);
  if (!pt) { p.hover.hidden = true; return; }
  // The value under the cursor in one dataset; null outside its grid.
  const read = (layer, dims) => {
    if (!layer.data) return null;
    const ex = dims[p.x].edges, ey = dims[p.y].edges, { values, counts, rows, cols } = layer.data;
    const col = Math.floor((pt[0] - ex[0]) / (ex[ex.length - 1] - ex[0]) * cols);
    const row = Math.floor((pt[1] - ey[0]) / (ey[ey.length - 1] - ey[0]) * rows);
    if (col < 0 || row < 0 || col >= cols || row >= rows) return null;
    const i = row * cols + col;
    return Number.isFinite(values[i]) ? `${fmtValue(values[i])} (${counts[i]} vox)` : 'no data';
  };
  const shown = compareShown();
  let text;
  if (!shown) {
    text = read(p, meta.dims);
    if (text === null) { p.hover.hidden = true; return; }
  } else {
    // Every dataset, the one under the cursor (or shown alone) first.
    const sets = [0, ...readyOthers().map((d) => d.slot + 1)];
    const [u0, u1, v0, v1] = p.view, um = (u0 + u1) / 2, vm = (v0 + v1) / 2;
    let under = shownIndex();
    if (shown === 'split') {
      under = others().length > 1 ? (pt[0] > um ? 1 : 0) + (pt[1] > vm ? 2 : 0)
        : (pt[0] - u0) / (u1 - u0) + (pt[1] - v0) / (v1 - v0) > 1 ? sets[1] : 0;
    }
    const order = sets.includes(under) ? [under, ...sets.filter((k) => k !== under)] : sets;
    text = order.map((k) => `${LETTERS[k]} ${read(layerOf(p, k), k ? slots[k - 1].meta.dims : meta.dims) ?? 'no data'}`).join(' · ');
  }
  const X = meta.dims[p.x], Y = meta.dims[p.y];
  p.hover.textContent = `${X.label} ${fmt(pt[0])}  ${Y.label} ${fmt(pt[1])}  →  ${text}`;
  p.hover.hidden = false;
}

function navigate(p, e) {
  const pt = pointAt(p, e);
  if (!pt) return;
  for (const q of panels) {
    const value = q.fixed === p.x ? pt[0] : q.fixed === p.y ? pt[1] : null;
    if (value === null) continue;
    q.center.value = q.slider.value = roundTo(clamp(value, q.lo, q.hi), q.step);
    paint(q.slider);
    request(q);
  }
}

function savePNG(p) {
  if (!p.data) return;
  const out = document.createElement('canvas'), s = uiScale();
  draw(p, out, Math.max(p.canvas.clientWidth / s, 480), Math.max(p.canvas.clientHeight / s, 360) + 24, 3, true);
  const X = meta.dims[p.x], Y = meta.dims[p.y], F = meta.dims[p.fixed];
  out.toBlob((blob) => download(blob, `${shownStem()}_${X.label}${Y.label}_${F.label}=${p.data.center}.png`));
}

/** File names of what the views show: one dataset, or all of them (see setStem). */
function shownStem() {
  const shown = compareShown();
  if (shown !== 'split') return shownIndex() ? stemOf(slots[shownIndex() - 1].name) : stem();
  return setStem([stem(), ...readyOthers().map((d) => stemOf(d.name))]);
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name.replace(/[^\w.=+-]+/g, '_');
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}

// ---- Symmetry --------------------------------------------------------------------------

const groupKey = (ops) => ops.map((m) => m.join()).sort().join(';');

function applySymmetry() {
  if (!panels.length) return; // viewer not ready yet
  let ops, maps;
  try {
    ops = closeGroup(parseOps($('sym-ops').value));
    maps = indexMaps(ops, meta.dims);
  } catch (err) {
    $('sym-status').textContent = err.message;
    $('sym-status').className = 'note error';
    return;
  }
  const key = groupKey(ops);
  const preset = PRESETS.find(([, gens]) => groupKey(closeGroup(parseOps(gens))) === key);
  $('sym-preset').value = preset ? preset[0] : 'custom';
  symmetry = { name: preset ? preset[0] : 'custom', ops, maps };
  for (const d of readyOthers()) Object.assign(d, compareMaps(d));
  if (readyOthers().length) showCompare();
  showSymmetry();
  describe();
  panels.forEach(request);
  requestIso();
  requestPowder();
}

function showSymmetry() {
  const { ops } = symmetry, statusEl = $('sym-status');
  statusEl.className = 'note';
  if (ops.length === 1) {
    statusEl.textContent = 'No symmetry averaging: each voxel is used as measured.';
  } else {
    let text = `${ops.length} operations; equivalent voxels are pooled with equal weight.`;
    const hkl = meta.dims.every((d) => d.basis) || meta.dims.every((d) => isDirect(d) && d.length > 0);
    if (meta.lattice && hkl) {
      const change = metricChange(ops, meta.lattice);
      if (change > 0.02) {
        text = `${ops.length} operations, but they change the cell metric by up to ${(100 * change).toFixed(0)}%. `
          + 'They are not symmetries of this lattice; check the setting.';
        statusEl.className = 'note warn';
      } else {
        text += ` Cell metric (${meta.lattice.source}) preserved to ${(100 * change).toFixed(2)}%.`;
      }
    } else if (!hkl) {
      text += ' Axes have no HKL basis or cell axes, so operations act on the display axes directly.';
    }
    statusEl.textContent = text;
  }
  $('sym-list').textContent = ops.map(formatOp).join('   ');
  $('sym-count').textContent = ops.length;
  updateStates();
}

// ---- Mask -------------------------------------------------------------------------------

function applyMask(clear = false) {
  if (!panels.length) return; // viewer not ready yet
  const radius = clear ? 0 : Math.round(Number($('mask-erode').value) || 0);
  const k = clear ? 0 : Number($('mask-k').value) || 0;
  if (radius < 0 || k < 0) {
    $('mask-status').className = 'note error';
    $('mask-status').textContent = 'Use a non-negative erosion radius and outlier cut.';
    return;
  }
  if (k > 0 && symmetry.ops.length < 3) {
    $('mask-status').className = 'note error';
    $('mask-status').textContent = 'The outlier cut compares symmetry equivalents: choose a Laue class with at least 3 operations first (e.g. 6/mmm).';
    return;
  }
  $('mask-apply').disabled = true;
  setJob('mask-a', 'Mask', radius || k ? 'Starting' : 'Clearing', 0);
  $('mask-state').textContent = 'working…';
  $('mask-state').className = 'state busy';
  $('mask-status').className = 'note';
  $('mask-status').textContent = radius || k ? 'Building mask…' : 'Clearing mask…';
  worker.postMessage({ type: 'mask', id: ++requestId, radius, k, maps: symmetry.maps, symmetry: symmetry.name });
  for (const d of readyOthers()) sendCompareMask(d, radius, k);
}

function showMask(seconds) {
  const statusEl = $('mask-status');
  statusEl.className = 'note';
  updateStates();
  if (!mask) {
    statusEl.textContent = 'No mask: all measured voxels are used.';
    return;
  }
  const { measured, edge, outlier, radius, k, group } = mask;
  const parts = [];
  if (radius) parts.push(`${pct(edge / measured)} within ${radius} voxel${radius === 1 ? '' : 's'} of coverage edges`);
  if (k) parts.push(`${pct(outlier / measured)} above ${k}σ of their ${group} equivalents`);
  statusEl.className = 'note ok';
  statusEl.textContent = `Removed ${pct((edge + outlier) / measured)} of measured voxels: ${parts.join(', ')}${seconds ? ` (${seconds.toFixed(1)} s)` : ''}.`;
}

// ---- Export for NEBULA3D ----------------------------------------------------------------

/** The Export card: what would be exported, or why it cannot be. */
function showExport() {
  if (exportJob || !meta) return;
  const statusEl = $('export-status');
  statusEl.title = '';
  let plan;
  try {
    plan = exportPlan(meta.dims, meta.lattice);
  } catch (err) {
    $('export-run').disabled = $('export-open').disabled = true;
    exportState('', 'unavailable');
    statusEl.className = 'note error';
    statusEl.textContent = err.message;
    return;
  }
  $('export-run').disabled = $('export-open').disabled = false;
  exportState('', 'ready');
  const sym = symmetry.ops.length > 1, voxels = plan.shape.reduce((a, b) => a * b);
  const parts = [
    `${plan.shape.join(' × ')} (H × K × L)${plan.padded ? ', padded to be symmetric about 0' : ''}`,
    sym ? `${symmetry.name} averaged` : 'not symmetrized',
    mask ? `mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : 'no mask',
  ];
  if (readyOthers().length) parts.push('dataset A');
  const warnings = [];
  if (!sym) warnings.push('NEBULA3D expects a symmetrized volume: choose a Laue class first.');
  if (voxels > 80e6) warnings.push(`NEBULA3D in the browser handles up to about 80 M voxels; this is ${Math.round(voxels / 1e6)} M, so use its desktop app.`);
  statusEl.className = warnings.length ? 'note warn' : 'note';
  statusEl.textContent = `${parts.join(' · ')}.${warnings.length ? ` ${warnings.join(' ')}` : ''}`;
}

function showExportProgress(label, fraction) {
  setJob('export', 'Export', label, fraction);
  $('export-status').className = 'note';
  $('export-status').textContent = `${label}… ${Math.round(100 * fraction)}%`;
}

/** Build the export in the worker: downloaded when done, or sent to NEBULA3D (`send`). */
function runExport(send = false) {
  if (!panels.length || exportJob) return false;
  let plan, maps;
  try {
    plan = exportPlan(meta.dims, meta.lattice);
    maps = symmetry.ops.length > 1 ? indexMaps(symmetry.ops, plan.paddedDims) : [IDENTITY_MAP];
  } catch (err) {
    $('export-status').className = 'note error';
    $('export-status').textContent = err.message;
    return false;
  }
  const sym = symmetry.ops.length > 1, { order, lo, size, shape, centers, ub, padded } = plan;
  const attrs = {
    source_file: sourceName,
    symmetry: sym ? symmetry.name : 'none',
    symmetry_ops: symmetry.ops.map(formatOp).join('; '),
    mask: mask ? `coverage-edge erosion ${mask.radius}, outlier cut ${mask.k} sigma, ${pct((mask.edge + mask.outlier) / mask.measured)} of measured voxels removed` : 'none',
    ub_source: meta.lattice.source,
    padded: padded ? 'yes' : 'no',
    created_by: 'NeXus Viewer, https://drthyang.github.io/neutron-nexus-viewer/',
    created: new Date().toISOString(),
  };
  exportJob = { name: `${stem()}_${sym ? `sym${symmetry.name.replace(/\//g, '')}` : 'unsym'}.nxs`, send, attrs };
  $('export-run').disabled = $('export-open').disabled = true;
  exportState('busy', 'working…');
  setJob('export', 'Export', 'Symmetrizing', 0);
  showExportProgress('Symmetrizing', 0);
  worker.postMessage({ type: 'export', id: ++requestId, plan: { order, lo, size, shape, centers, ub }, maps, attrs });
  return true;
}

/** Show the export's state and a message; `file`, the full name a message shortens, is its tooltip. */
function exportNote(kind, state, text, file = '') {
  exportState(kind, state);
  $('export-status').className = `note${kind === 'ok' ? ' ok' : kind === 'error' ? ' error' : ''}`;
  $('export-status').textContent = text;
  $('export-status').title = file;
}

/** The export's state badge, its section summary and the pipeline's NEBULA3D step. */
function exportState(kind, state) {
  $('export-state').textContent = state;
  $('export-state').className = `state ${kind}`.trim();
  $('sum-export').textContent = `I(Q) · NEBULA3D · ${state}`;
  $('pipe-export').className = kind === 'ok' ? 'ok' : 'off';
  $('pipe-export-text').textContent = { sent: 'volume sent to NEBULA3D', saved: 'volume saved for NEBULA3D' }[state] ?? 'next analysis: export the volume';
}

/**
 * Open in NEBULA3D: open its page in a new tab (now, while this is a click),
 * build the volume meanwhile, and post the file once NEBULA3D says it is
 * ready. The protocol is described in NEBULA3D's web/src/api/importHandoff.ts.
 *
 * On this origin (the deployed pair) the tab opens with `noopener` and the
 * messages go over a BroadcastChannel: tabs of one site that hold a window
 * reference to each other share a renderer process, so reloading or closing
 * this page would also stop a run in NEBULA3D. Another origin (a local dev
 * server) needs the window reference for postMessage.
 */
function openInNebula() {
  if (!panels.length || exportJob) return;
  const url = new URL(NEBULA3D_URL, location.href), id = crypto.randomUUID();
  url.searchParams.set('import', 'nexus-viewer');
  url.searchParams.set('id', id);
  url.searchParams.set('from', location.origin);
  const separate = url.origin === location.origin && typeof BroadcastChannel === 'function';
  // With noopener, window.open returns null even when the tab opens; watchHandoff reports a blocked one.
  const win = window.open(url.href, '_blank', separate ? 'noopener' : '');
  if (!separate && !win) {
    exportNote('error', 'blocked', 'The browser blocked the new tab: allow pop-ups for this site, or use Download.');
    return;
  }
  endHandoff();
  const channel = separate ? new BroadcastChannel(`nebula3d-import:${id}`) : null;
  channel?.addEventListener('message', (e) => onHandoffMessage(e.data));
  handoff = { id, win, channel, origin: url.origin, ready: false, file: null, meta: null, sent: false, since: Date.now() };
  handoff.watch = setInterval(watchHandoff, 1000);
  if (!runExport(true)) endHandoff('the volume could not be built.');
}

/** Post `message` to the NEBULA3D tab of handoff `h`. */
function postHandoff(h, message) {
  if (h.channel) h.channel.postMessage(message);
  else h.win.postMessage(message, h.origin);
}

/** Show the waiting NEBULA3D tab how far the volume is built (once it is listening). */
function forwardProgress(label, fraction) {
  const h = handoff;
  if (!h || h.file) return;
  h.progress = { label, fraction };
  if (!h.ready) return;
  try { postHandoff(h, { type: 'nebula3d-import-progress', id: h.id, label, fraction }); } catch { /* tab gone */ }
}

function sendHandoff() {
  const h = handoff;
  if (!h?.ready || !h.file || h.sent) return;
  h.sent = true;
  postHandoff(h, { type: 'nebula3d-import', id: h.id, schema: 'nexus-viewer/1', file: h.file, meta: h.meta });
  exportNote('busy', 'sending…', `Sent ${shortName(h.file.name)} to NEBULA3D; loading it there…`, h.file.name);
}

/** Stop the handoff; with a `reason`, tell NEBULA3D (so it stops waiting) and show it. */
function endHandoff(reason = null) {
  const h = handoff;
  handoff = null;
  if (!h) return;
  clearInterval(h.watch);
  if (reason) {
    try { postHandoff(h, { type: 'nebula3d-import-cancel', id: h.id, message: reason }); } catch { /* tab gone */ }
    exportNote('error', 'not sent', `Not sent to NEBULA3D: ${reason}`);
  }
  h.channel?.close();
}

function watchHandoff() {
  const h = handoff;
  if (!h) return;
  if (h.win?.closed) endHandoff('the NEBULA3D tab was closed.');
  else if (h.channel && !h.ready && h.file && Date.now() - h.since > 20_000) {
    // NEBULA3D answers as soon as its page loads, so this is most likely a blocked tab.
    exportNote('error', 'no answer', 'NEBULA3D has not answered. If no tab opened, allow pop-ups for this site and try again, or use Download and Load volume… there.');
  } else if (!h.ready && h.file && Date.now() - h.since > 120_000) {
    exportNote('error', 'no answer', 'NEBULA3D has not answered after 2 minutes. Check its tab, or use Download and Load volume… there.');
  }
}

/** A message from the NEBULA3D tab, by window or channel. */
function onHandoffMessage(data) {
  const h = handoff;
  if (!h || data?.id !== h.id) return;
  if (data.type === 'nebula3d-import-ready') {
    // On the first "ready", show NEBULA3D the progress it missed while loading.
    const first = !h.ready;
    h.ready = true;
    if (first && h.progress) forwardProgress(h.progress.label, h.progress.fraction);
    sendHandoff();
  } else if (data.type === 'nebula3d-import-loaded') {
    const name = h.file?.name ?? 'the volume';
    endHandoff();
    exportNote('ok', 'sent', `NEBULA3D loaded ${shortName(name)} and selected it as its dataset; continue in its tab.`, name);
  } else if (data.type === 'nebula3d-import-error') {
    endHandoff();
    exportNote('error', 'failed', `NEBULA3D could not load the volume: ${data.message}`);
  }
}

addEventListener('message', (e) => {
  const h = handoff;
  if (h?.win && e.source === h.win && e.origin === h.origin) onHandoffMessage(e.data);
});

// ---- 3-D view ---------------------------------------------------------------------------

async function setup3D() {
  const current = iso;
  try {
    const { View3D } = await import('./view3d.js');
    if (iso !== current) return; // another file was opened meanwhile
    view3d?.renderer.dispose();
    view3d = new View3D(iso.canvas);
  } catch (err) {
    iso.caption.textContent = `3-D view unavailable: ${err.message}`;
    return;
  }
  // An empty level returns to the automatic (99.5th percentile) level.
  iso.levelInput.onchange = () => {
    iso.userLevel = iso.levelInput.value.trim() !== '';
    requestIso();
  };
  iso.slider.oninput = () => {
    if (!iso.range) return;
    iso.userLevel = true;
    iso.levelInput.value = sig(sliderToLevel(Number(iso.slider.value)), 3);
    requestIso();
  };
  iso.opacity.oninput = () => view3d.setOpacity(Number(iso.opacity.value));
  iso.sliceOpacity.oninput = () => view3d.setSliceOpacity(Number(iso.sliceOpacity.value));
  iso.reset.onclick = () => view3d.resetView();
  iso.save.onclick = () => view3d.snapshot((blob) => download(blob, `${stem()}_3d.png`));
  update3D();
  requestIso();
}

// The slider spans the median to the 99.99th percentile of positive block means, logarithmically.
function sliderToLevel(t) {
  const [lo, hi] = iso.range;
  return lo * (hi / lo) ** (t / 1000);
}

function levelToSlider(level) {
  const [lo, hi] = iso.range;
  return clamp(1000 * Math.log(level / lo) / Math.log(hi / lo), 0, 1000);
}

function requestIso() {
  if (!iso || !view3d) return;
  const level = iso.userLevel ? Number(iso.levelInput.value) : null;
  if (level !== null && !Number.isFinite(level)) return;
  iso.wanted = { maxBins: Number($('iso-grid').dataset.value), level, ops: symmetry.ops, symmetry: symmetry.name };
  if (!iso.busy) sendIso();
}

function sendIso() {
  const q = iso.wanted;
  iso.wanted = null;
  iso.busy = true;
  worker.postMessage({ type: 'iso', id: ++requestId, ...q });
}

/** Display -> Cartesian reciprocal space (see cartesianBasis()), with the cell angles drawn. */
function displayBasis() {
  const cell = meta.lattice && (settings.angles === 'nominal' ? nominalCell(meta.lattice) : meta.lattice);
  return cartesianBasis(meta.dims, cell);
}

function update3D() {
  if (!view3d || !settings) return;
  const { T } = displayBasis();
  view3d.setFrame(T, meta.dims.map((d) => viewRange(d, settings.limit)), meta.dims.map((d) => d.label));
  const slices = panels.filter((p) => p.data).map((p) => {
    const X = meta.dims[p.x], Y = meta.dims[p.y], u = viewRange(X, settings.limit), v = viewRange(Y, settings.limit);
    return {
      fixed: p.fixed, x: p.x, y: p.y, center: p.data.center, ...sliceTexture(p, settings, u, v), u, v,
      ex: [X.edges[0], X.edges[X.edges.length - 1]], ey: [Y.edges[0], Y.edges[Y.edges.length - 1]],
    };
  });
  view3d.setSlices(slices, $('iso-slices').checked);
  show3DCut();
}

/** The line cut in the 3-D view: its line through the slices, inside the rod of voxels it averages. */
function show3DCut() {
  if (!view3d) return;
  let spec = null;
  try {
    if (cut?.line) spec = cutSpec();
  } catch { /* the cut view says why */ }
  if (!spec) {
    view3d.setCut(null);
    return;
  }
  const { a, b, dom, edges, radius } = spec, along = (s) => a.map((ad, d) => ad + ((s - a[dom]) / (b[dom] - a[dom])) * (b[d] - ad));
  // Where the cut passes through the slices it does not lie in.
  const pierce = panels.filter((q) => q.data).flatMap((q) => {
    const d = q.fixed, t = (q.data.center - a[d]) / (b[d] - a[d]);
    return t > 0 && t < 1 ? [{ at: a.map((ad, i) => ad + t * (b[i] - ad)), axis: d }] : [];
  });
  view3d.setCut({ line: [a, b], rod: [along(edges[0]), along(edges[edges.length - 1])], radius, pierce });
}

/**
 * The texture of a slice plane in the 3-D view, on A's pixel grid. When
 * comparing, it matches the 2-D view: the datasets' parts of the visible range
 * (u, v), see splitLayout(), or one dataset alone.
 */
function sliceTexture(p, s, u, v) {
  const shown = compareShown();
  if (!shown || shown === 'a') return { image: layerImage(p, s), key: p.imageKey };
  const split = splitLayout(u[0], u[1], v[0], v[1]);
  const images = split.parts.map(([k, part]) => [k, part, layerOf(p, k).data ? layerImage(layerOf(p, k), s) : null]);
  const key = `${images.map(([k, , image]) => (image ? layerOf(p, k).imageKey : '-')).join('|')}|${shown}|${others().length}|${u}|${v}`;
  if (p.textureKey === key) return { image: p.texture, key };
  const { rows, cols } = p.data, ex = meta.dims[p.x].edges, ey = meta.dims[p.y].edges;
  const dx = (ex[ex.length - 1] - ex[0]) / cols, dy = (ey[ey.length - 1] - ey[0]) / rows;
  const px = (uu, vv) => [(uu - ex[0]) / dx, (vv - ey[0]) / dy];
  const canvas = p.texture ?? document.createElement('canvas');
  canvas.width = cols;
  canvas.height = rows;
  const c = canvas.getContext('2d');
  const clip = (pts) => {
    c.beginPath();
    pts.forEach(([x, y], i) => (i ? c.lineTo(x, y) : c.moveTo(x, y)));
    c.closePath();
    c.clip();
  };
  c.imageSmoothingEnabled = false;
  for (const [k, part, image] of images) {
    if (!image) continue;
    // Each dataset's image placed on A's pixel grid, clipped to its part.
    const dims = k ? slots[k - 1].meta.dims : meta.dims, bx = dims[p.x].edges, by = dims[p.y].edges, { rows: rb, cols: cb } = layerOf(p, k).data;
    c.save();
    if (shown === 'split') clip(part.map(([uu, vv]) => px(uu, vv)));
    c.transform((bx[bx.length - 1] - bx[0]) / cb / dx, 0, 0, (by[by.length - 1] - by[0]) / rb / dy, (bx[0] - ex[0]) / dx, (by[0] - ey[0]) / dy);
    c.imageSmoothingEnabled = false;
    c.drawImage(image, 0, 0);
    c.restore();
  }
  if (split.lines.length) {
    c.strokeStyle = '#ffffff';
    c.lineWidth = Math.max(1, Math.max(cols, rows) / 120);
    c.beginPath();
    for (const [from, to] of split.lines) { c.moveTo(...px(...from)); c.lineTo(...px(...to)); }
    c.stroke();
  }
  p.texture = canvas;
  p.textureKey = key;
  return { image: canvas, key };
}

/** Apply the symmetry (a Laue class name) and mask ("r" or "r,k") requested for this file. */
function applyPendingProcessing() {
  const p = pendingProcessing;
  pendingProcessing = null;
  const preset = p?.sym && PRESETS.find(([name]) => name === p.sym);
  if (preset) {
    $('sym-preset').value = preset[0];
    $('sym-ops').value = preset[0] === '1' ? '' : preset[1];
    applySymmetry();
  }
  if (p?.mask) {
    const [radius = 0, k = 0] = p.mask.split(',').map(Number);
    $('mask-erode').value = radius;
    $('mask-k').value = k;
    applyMask();
  }
}

/** The example: a synthetic crystal at 300 K (A) and 10 K (B), symmetrized and masked. */
function openDemo() {
  pendingCompare = [DEMO.compare];
  pendingProcessing = { sym: DEMO.sym, mask: DEMO.mask };
  openURL(DEMO.url);
}

// ---- I(Q) ---------------------------------------------------------------------------------

// Curve colors of datasets A to D: the accent, the H- and L-axis hues, and violet.
const CURVE = { A: '#2f74e6', B: '#d98a0b', C: '#13a36b', D: '#8e44c9' };
const newCurveLayer = () => ({ data: null, key: '', busy: false, wanted: null, error: '', progress: null });

/**
 * Recompute I(Q): now if a download waits, or if it was asked for (Show I(Q), the
 * I(Q) switch or Compute) and its view is shown; otherwise once that holds.
 */
function requestPowder() {
  if (!powder) return;
  powder.stale = true;
  updatePowder();
}

/**
 * Send the I(Q) requests when they are due. The other datasets get A's shells,
 * which by default reach the farthest of the grids; their masks are built first.
 */
function updatePowder() {
  if (!powder?.stale || !((powder.requested && plotShown(powder)) || powder.download)) return;
  powder.stale = false;
  const read = (input, zero = false) => {
    const text = input.value.trim(), x = Number(text);
    if (text === '') return null;
    if (!(x > 0 || (zero && x === 0))) throw new Error(`${input.getAttribute('aria-label')} must be a positive number, or empty for automatic.`);
    return x;
  };
  const split = Number($('iq-split').dataset.value);
  // Per other dataset: its |Q| extent, then its plan on A's shells, or why it has none.
  const more = readyOthers().map((d) => ({ d, extent: null, plan: null, error: '' }));
  let plan, shells;
  try {
    const bins = parseBins(powder.bins.value), qmin = read(powder.qmin, true), qmax = read(powder.qmax);
    let top = qExtent(meta.dims, meta.lattice).top;
    for (const o of more) {
      try {
        o.extent = qExtent(o.d.meta.dims, o.d.meta.lattice);
        top = Math.max(top, o.extent.top);
      } catch (err) {
        o.error = err.message;
      }
    }
    plan = powderPlan(meta.dims, meta.lattice, { bins, qmin, qmax: qmax ?? top, split });
    shells = JSON.stringify([plan.bins, qmin, qmax ?? top, split]);
    for (const o of more) if (o.extent) o.plan = powderPlan(o.d.meta.dims, o.d.meta.lattice, { edges: plan.edges, split });
  } catch (err) {
    Object.assign(powder.a, { data: null, error: err.message, key: '' });
    powder.plan = null;
    powder.download = false;
    showPowder();
    drawPlot(powder);
    return;
  }
  powder.plan = plan;
  const group = groupKey(symmetry.ops);
  queuePowder(powder.a, worker, { plan, maps: symmetry.maps, symmetry: symmetry.name }, `${shells}|${group}|${maskVersion}`);
  for (const { d, plan: own, error: why } of more) {
    if (own && !d.maskBusy) {
      queuePowder(powder.more[d.slot], d.worker, { plan: own, maps: d.maps, symmetry: d.maps.length > 1 ? symmetry.name : '1' },
        `${shells}|${group}|${d.maps.length}|${d.maskVersion ?? 0}`);
    } else if (why) {
      Object.assign(powder.more[d.slot], { data: null, error: why, key: '' });
    }
  }
  showPowder();
  finishPowderDownload();
}

// The ΔQ slider sets a single step on a log scale: a width from a tenth of the
// shortest bin step in |Q| to 20 times it, or ΔQ/Q from 0.1% to 20% when the
// step is negative (logarithmic bins). Ranges with their own steps are typed.
const binsRange = (log) => (log ? [0.001, 0.2] : [powder.step / 10, powder.step * 20]);

/** The Rebin parameters in the ΔQ field, or the automatic step; null while the field cannot be read. */
function currentBins() {
  try {
    return parseBins(powder.bins.value) ?? [Number(powder.step.toPrecision(1))];
  } catch {
    return null;
  }
}

/** Place the slider at the typed step; with ranges, Q min and Q max come from the field and the slider rests. */
function syncBinsSlider() {
  const bins = powder.step ? currentBins() : null, single = bins?.length === 1;
  powder.slider.disabled = !single;
  powder.qmin.disabled = powder.qmax.disabled = bins?.length > 1;
  powder.qmin.placeholder = single && bins[0] < 0 ? String(Number(powder.step.toPrecision(2))) : '0';
  if (!single) {
    powder.slider.title = bins ? 'These bins have ranges with their own steps: edit them in the field' : '';
    return;
  }
  const log = bins[0] < 0, [lo, hi] = binsRange(log);
  powder.slider.value = 1000 * clamp(Math.log(Math.abs(bins[0]) / lo) / Math.log(hi / lo), 0, 1);
  powder.slider.title = `${binsLabel(bins)}: drag to change it, recomputing as you go`;
  paint(powder.slider);
}

function slideBins() {
  const bins = currentBins(), log = bins?.length === 1 && bins[0] < 0, [lo, hi] = binsRange(log);
  const step = sig(lo * (hi / lo) ** (Number(powder.slider.value) / 1000), 2);
  powder.bins.value = String(log ? -step : step);
  powder.slider.title = `${binsLabel([log ? -step : step])}: drag to change it, recomputing as you go`;
  requestPowder();
}

/** Q bins in words: "ΔQ 0.05 Å⁻¹", "ΔQ/Q 1%", or the Rebin parameters of ranges. */
function binsLabel(bins) {
  if (bins.length > 1) return `bins ${bins.join(', ')}`;
  return bins[0] > 0 ? `ΔQ ${fmt(bins[0], 4)} Å⁻¹` : `ΔQ/Q ${sig(-100 * bins[0], 2)}%`;
}

/** The last shell with data, and the |Q| where it ends. */
function lastShell(data) {
  let b = data.voxels.length - 1;
  while (b > 0 && !(data.voxels[b] > 0)) b--;
  return b;
}
const powderEnd = (data) => data.edges[lastShell(data) + 1];
/** The centre of shell b (the midpoint of its edges). */
const shellMid = (data, b) => (data.edges[b] + data.edges[b + 1]) / 2;

/** The shell of `edges` that holds q, or -1. */
function shellAt(edges, q) {
  if (!(q >= edges[0] && q < edges[edges.length - 1])) return -1;
  let b = 0, e = edges.length - 1;
  while (e - b > 1) {
    const m = (b + e) >> 1;
    if (q >= edges[m]) b = m;
    else e = m;
  }
  return b;
}

/** Ask a dataset's worker for I(Q), unless it has (or is computing) the same `key`. One request in flight per dataset. */
function queuePowder(layer, w, request, key) {
  if (key === layer.key && !layer.error) return;
  layer.key = key;
  layer.wanted = request;
  if (!layer.busy) sendPowder(layer, w);
}

function sendPowder(layer, w) {
  const q = layer.wanted, k = layer === powder.a ? 0 : powder.more.indexOf(layer) + 1;
  layer.wanted = null;
  layer.busy = true;
  layer.progress = null;
  setJob(`iq-${'abcd'[k]}`, k ? `I(Q) ${LETTERS[k]}` : 'I(Q)', 'Starting', 0);
  w.postMessage({ type: 'powder', id: ++requestId, ...q });
}

/** I(Q) of dataset k (0 for A) has arrived, or failed with `message`. */
function powderResult(k, msg, message = '') {
  if (!powder) return;
  const layer = layerOf(powder, k);
  layer.busy = false;
  layer.progress = null;
  Object.assign(layer, msg ? { data: msg, error: '' } : { data: null, error: message });
  if (layer.wanted) sendPowder(layer, workerOf(k));
  else endJob(`iq-${'abcd'[k]}`);
  showPowder();
  drawPlot(powder);
  finishPowderDownload();
}

function powderProgress(k, { label, fraction }) {
  const layer = powder && layerOf(powder, k);
  if (!layer?.busy) return;
  layer.progress = { label, fraction };
  setJob(`iq-${'abcd'[k]}`, k ? `I(Q) ${LETTERS[k]}` : 'I(Q)', label, fraction);
  showPowder();
}

/** The I(Q) card in the Export section, and the view's header. */
function showPowder() {
  if (!powder) return;
  const statusEl = $('powder-status');
  const state = (kind, text) => {
    $('powder-state').textContent = text;
    $('powder-state').className = `state ${kind}`.trim();
  };
  const layers = [['A', powder.a], ...readyOthers().map((d) => [d.letter, powder.more[d.slot]])];
  const working = layers.filter(([, l]) => l.busy);
  powder.compute.hidden = powder.requested || !!powder.a.error;
  if (working.length) {
    const fraction = working.reduce((s, [, l]) => s + (l.progress?.fraction ?? 0), 0) / working.length;
    const label = working.find(([, l]) => l.progress)?.[1].progress.label ?? 'Starting';
    state('busy', 'working…');
    statusEl.className = 'note';
    statusEl.textContent = `${label}${layers.length > 1 ? ` (${working.map(([k]) => k).join(', ')})` : ''}… ${Math.round(100 * fraction)}%`;
    powder.caption.textContent = `computing… ${Math.round(100 * fraction)}%`;
    powder.caption.title = '';
    return;
  }
  const A = powder.a.data;
  if (!A) {
    state('', powder.a.error ? 'unavailable' : 'off');
    statusEl.className = powder.a.error ? 'note error' : 'note';
    statusEl.textContent = powder.a.error
      || (powder.requested ? 'Computed when its view is shown.' : 'Not computed yet: choose Show I(Q), or I(Q) in the header of the 3-D view.');
    powder.caption.textContent = powder.a.error ? 'unavailable' : '';
    powder.caption.title = powder.a.error;
    return;
  }
  const sets = layers.filter(([, l]) => l.data);
  const last = Math.max(...sets.map(([, l]) => lastShell(l.data))), end = A.edges[last + 1];
  const errors = sets.map(([k, l]) => [k, l.data.errors ? `σ from ${l.data.errors}` : l.data.errorsNote || 'no uncertainties in the file']);
  const facts = [
    `${last + 1} shells (${binsLabel(powder.plan.bins)}) from ${fmt(A.edges[0], 3)} to ${fmt(end, 2)} Å⁻¹`,
    A.order > 1 ? `${A.symmetry} averaged` : 'not symmetrized',
    A.masked && mask ? `mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : 'no mask',
    ...(errors.every(([, t]) => t === errors[0][1]) ? [errors[0][1].replace('the file', sets.length > 1 ? 'the files' : 'the file')] : errors.map(([k, t]) => `${k}: ${t}`)),
  ];
  const warnings = layers.slice(1).filter(([, l]) => l.error).map(([k, l]) => `${k}: ${l.error}`);
  if (powder.stale) warnings.push('Settings changed: it is recomputed when its view is shown.');
  state(powder.stale ? '' : 'ok', powder.stale ? 'out of date' : 'ready');
  statusEl.className = warnings.length ? 'note warn' : 'note';
  statusEl.textContent = `${sets.length > 1 ? `${listText(sets.map(([k]) => k))}: ` : ''}${facts.join(' · ')}.${warnings.length ? ` ${warnings.join(' ')}` : ''}`;
  powder.caption.textContent = `${A.split > 1 ? `${A.split}³ split` : 'centres'}${A.order > 1 ? ` · ${A.symmetry}` : ''}${A.masked ? ' · masked' : ''}${A.errors ? ' · ±σ' : ''}`;
  powder.caption.title = `${facts.join('\n')}\nVoxels split into ${A.split}³ sub-cells · ${sets.map(([k, l]) => `${sets.length > 1 ? `${k} ` : ''}${l.data.seconds.toFixed(1)} s`).join(', ')}`;
}

/** Show the I(Q) view: in the fourth place, or as the single view. */
function showPowderView() {
  if (!powder) return;
  powder.requested = true;
  if (layout === 'single') setLayout('single', 'iq');
  else setSlot('iq');
  views.iq.section.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
}

/** Save I(Q) as text, computing it first if needed. */
function downloadPowder() {
  if (!powder) return;
  powder.requested = powder.download = true;
  updatePowder();
  finishPowderDownload();
}

function finishPowderDownload() {
  if (!powder?.download || powder.stale || others().some((d) => d.maskBusy)) return;
  if ([powder.a, ...readyOthers().map((d) => powder.more[d.slot])].some((l) => l.busy || l.wanted)) return;
  powder.download = false;
  if (!powder.a.data) return;
  download(new Blob([powderText()], { type: 'text/plain' }), `${setStem(curveSets(powder).map(([, name]) => stemOf(name)))}_IQ.dat`);
}

/** The datasets with a curve in `plot` (I(Q) or the cut), A first: [letter, file name, data]. */
function curveSets(plot) {
  return [['A', sourceName, plot.a.data], ...readyOthers().map((d) => [d.letter, d.name, plot.more[d.slot].data])].filter(([, , data]) => data);
}

/** I(Q) as text: a commented header, then Q and, per dataset, I, σ, coverage and voxels. */
function powderText() {
  const sets = curveSets(powder);
  const both = sets.length > 1, { edges, split, frame } = sets[0][2], l = meta.lattice, bins = powder.plan.bins;
  const rows = Math.max(...sets.map(([, , d]) => lastShell(d))) + 1;
  // The same shells as Rebin parameters, with a single step's range written out.
  const rebin = bins.length > 1 ? bins : [edges[0], bins[0], edges[edges.length - 1]];
  const num = (x) => (Number.isFinite(x) ? String(Number(x.toPrecision(7))) : 'nan');
  const columns = both ? sets.flatMap(([k]) => [`I_${k}`, `sigma_${k}`, `coverage_${k}`, `voxels_${k}`]) : ['I', 'sigma', 'coverage', 'voxels'];
  const lines = [
    `# I(Q), the powder average of ${listText(sets.map(([k, name]) => (both ? `${k} = ${name}` : name)))}`,
    `# Written by NeXus Viewer (https://drthyang.github.io/neutron-nexus-viewer/) on ${new Date().toISOString()}`,
    frame === 'HKL'
      ? `# |Q| in 1/Angstrom, with 2*pi, from the cell (${l.source}): a b c = ${[l.a, l.b, l.c].map((x) => x.toFixed(4)).join(' ')} Angstrom, alpha beta gamma = ${[l.alpha, l.beta, l.gamma].map((x) => x.toFixed(3)).join(' ')} deg`
      : '# |Q| in 1/Angstrom, from the Q axes',
    `# Shells: ${binsLabel(bins).replace('Å⁻¹', '1/Angstrom')}; as Mantid Rebin parameters ${rebin.map((x) => Number(x.toPrecision(10))).join(', ')}`,
    `#   (a negative step is logarithmic; a range ends with a bin of 0.25 to 1.25 steps). Q is the shell centre, the midpoint of its edges.`,
    `# Voxels split into ${split}^3 sub-cells, each binned by its own |Q|`,
    `# Symmetry: ${symmetry.ops.length > 1 ? `${symmetry.name} (${symmetry.ops.length} operations), equivalent voxels pooled` : 'none'}`,
    `# Mask: ${mask ? `coverage-edge erosion ${mask.radius}${mask.k ? `, outlier cut ${mask.k} sigma` : ''}${both ? ', built for each dataset' : `, ${pct((mask.edge + mask.outlier) / mask.measured)} of measured voxels removed`}` : 'none'}`,
    '# I: the intensity integrated over the part of the shell with data, divided by the volume of that part',
    '#    (unmeasured and masked voxels are left out, not counted as zero; each symmetry orbit counts with its',
    '#    multiplicity). I * coverage is the same integral divided by the volume of the whole shell.',
    `# sigma: ${sets.map(([k, , d]) => `${both ? `${k} ` : ''}${d.errors ? `propagated from ${d.errors}` : 'nan, no uncertainties in the file'}`).join('; ')}`,
    '# coverage: fraction of the shell volume with data; voxels: number of voxels with data in the shell',
    `# Q ${columns.join(' ')}`,
  ];
  for (let b = 0; b < rows; b++) {
    const cells = [((edges[b] + edges[b + 1]) / 2).toFixed(6)];
    for (const [, , d] of sets) cells.push(num(d.intensity[b]), num(d.sigma[b]), d.coverage[b].toFixed(5), num(d.voxels[b]));
    lines.push(cells.join(' '));
  }
  return `${lines.join('\n')}\n`;
}

// ---- 1-D plots: I(Q) and the line cut -----------------------------------------------------
//
// Both are drawn, read and zoomed by the functions below from a plot object:
// its canvas, hover chip and Reset zoom button, its `a` and `b` layers (one per
// dataset, each with `data` = { edges, intensity, sigma }), the window of the
// last drawing (`view`), the zoom, a drag and the bin under the cursor (`at`), and
//   log(), band(), cover()   the intensity scale, the ±σ band and the coverage
//                            line on a right-hand axis (I(Q) only)
//   end(data)                where the x range ends: the last bin with data
//   xMin                     the smallest x a click zoom may start at
//   xLabel(), yLabel(log)    the axis titles
//   empty()                  what a plot without curves says, { text, error }, or null
//   title, extras(data)      an export's title and its line of settings
//   readout(b, layers)       the hover text for bin b
//   onHover()                optional: called when `at` changes

/** The curves of a plot shown: every dataset when comparing, or the one shown alone (following A / Split / B …). */
function shownLayers(plot) {
  const shown = compareShown(), only = shown && shown !== 'split' ? shownIndex() : null;
  return curveSets(plot).filter(([letter]) => only === null || LETTERS[only] === letter)
    .map(([letter, name, data]) => ({ letter, name, data, color: CURVE[letter] }));
}

/**
 * The plot window: the zoom, or all bins with data and the range of the
 * curves (and bands) over them. Intensities are in plot units: log10(I) on a
 * log scale.
 */
function curveWindow(plot, layers, log, band) {
  const [x0, x1] = plot.zoom?.x ?? [Math.min(...layers.map((l) => l.data.edges[0])), Math.max(...layers.map((l) => plot.end(l.data)))];
  if (plot.zoom?.y) return { x0, x1, y0: plot.zoom.y[0], y1: plot.zoom.y[1] };
  let lo = Infinity, hi = -Infinity;
  for (const { data } of layers) {
    for (let b = 0; b < data.intensity.length; b++) {
      const q = shellMid(data, b), v = data.intensity[b], s = band && Number.isFinite(data.sigma[b]) ? data.sigma[b] : 0;
      if (q < x0 || q > x1 || !Number.isFinite(v) || (log && v <= 0)) continue;
      lo = Math.min(lo, log && v - s <= 0 ? v : v - s);
      hi = Math.max(hi, v + s);
    }
  }
  if (!(hi >= lo)) [lo, hi] = log ? [1, 10] : [0, 1];
  if (log) {
    const a = Math.log10(lo), b = Math.log10(hi), pad = Math.max(0.05 * (b - a), 0.1);
    return { x0, x1, y0: a - pad, y1: b + pad };
  }
  // Intensities start from 0 unless some are negative.
  const base = Math.min(0, lo), span = hi - base || Math.abs(hi) || 1;
  return { x0, x1, y0: base - (lo < 0 ? 0.05 * span : 0), y1: hi + 0.06 * span };
}

/**
 * Ticks of a log axis from 10^a to 10^b, as exponents: decades, with 2 and 5
 * on short spans, and evenly spaced values on spans shorter than that.
 */
function logTicks(a, b) {
  const out = [];
  for (let k = Math.floor(a); k <= Math.ceil(b); k++) {
    for (const m of b - a < 2.5 ? [1, 2, 5] : [1]) {
      const t = k + Math.log10(m);
      if (t >= a && t <= b) out.push(t);
    }
  }
  if (out.length < 3) return niceTicks(10 ** a, 10 ** b, 4).filter((v) => v > 0).map(Math.log10);
  const step = Math.ceil(out.length / 7);
  return out.filter((_, i) => i % step === 0);
}

/** Fill `text` centred at (x, y), wrapped to `width`. */
function wrapText(c, text, x, y, width, lineHeight) {
  const lines = [];
  for (const word of text.split(' ')) {
    const longer = lines.length ? `${lines[lines.length - 1]} ${word}` : word;
    if (lines.length && c.measureText(longer).width <= width) lines[lines.length - 1] = longer;
    else lines.push(word);
  }
  lines.forEach((line, i) => c.fillText(line, x, y + (i - (lines.length - 1) / 2) * lineHeight));
}

/** Whether a plot's view is on screen. */
const plotShown = (plot) => !!plot && plot.canvas.clientWidth > 0;

function drawPlot(plot) {
  if (!plotShown(plot)) return;
  const s = uiScale();
  drawCurves(plot, plot.canvas, plot.canvas.clientWidth / s, plot.canvas.clientHeight / s, s * (devicePixelRatio || 1));
}

/** The pointer on a plot's canvas, in its drawing coordinates. */
function plotLocal(plot, e) {
  const r = plot.canvas.getBoundingClientRect(), s = uiScale();
  return [(e.clientX - r.left) / s, (e.clientY - r.top) / s];
}

/**
 * Draw a plot into `canvas` (w × h CSS px): curves with ±σ bands, and for I(Q)
 * the shell coverage on the right-hand axis. Exports add a title and a legend.
 */
function drawCurves(plot, canvas, w, h, dpr, exporting = false) {
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const c = canvas.getContext('2d');
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  c.fillStyle = '#ffffff';
  c.fillRect(0, 0, w, h);
  if (w < 80 || h < 60) return;
  const layers = shownLayers(plot);
  if (!layers.length) {
    if (!exporting) plot.view = null;
    const note = plot.empty();
    if (!note) return;
    c.fillStyle = note.error ? '#d64545' : INK2;
    c.font = `12.5px ${SANS}`;
    c.textAlign = 'center';
    c.textBaseline = 'middle';
    wrapText(c, note.text, w / 2, h / 2, Math.min(w - 40, 420), 18);
    return;
  }
  const log = plot.log(), band = plot.band(), cover = plot.cover();
  const { x0, x1, y0, y1 } = curveWindow(plot, layers, log, band);
  const pad = { l: 62, r: cover ? 52 : 18, t: exporting ? 46 : 14, b: 44 };
  const aw = Math.max(10, w - pad.l - pad.r), ah = Math.max(10, h - pad.t - pad.b), bottom = pad.t + ah;
  const X = (q) => pad.l + ((q - x0) / (x1 - x0)) * aw;
  const T = (t) => bottom - ((t - y0) / (y1 - y0)) * ah;
  const Y = (v) => T(log ? Math.log10(v) : v);
  // Coverage: 0 at the bottom, 100% a little below the top.
  const C = (f) => bottom - (f / 1.08) * ah;
  if (!exporting) plot.view = { x0, x1, y0, y1, pad, aw, ah };
  const xt = niceTicks(x0, x1, Math.max(2, Math.round(aw / 80)));
  const yt = log ? logTicks(y0, y1) : niceTicks(y0, y1, Math.max(2, Math.round(ah / 55)));
  c.strokeStyle = '#eef1f4';
  c.lineWidth = 1;
  c.beginPath();
  for (const t of xt) { const x = Math.round(X(t)) + 0.5; c.moveTo(x, pad.t); c.lineTo(x, bottom); }
  for (const t of yt) { const y = Math.round(T(t)) + 0.5; c.moveTo(pad.l, y); c.lineTo(pad.l + aw, y); }
  c.stroke();

  c.save();
  c.beginPath();
  c.rect(pad.l, pad.t, aw, ah);
  c.clip();
  // Bins with data, in runs: a bin without data breaks the curve.
  const runs = (data) => {
    const out = [];
    let run = null;
    data.intensity.forEach((v, b) => {
      if (Number.isFinite(v) && (!log || v > 0)) (run ??= []).push(b);
      else if (run) { out.push(run); run = null; }
    });
    if (run) out.push(run);
    return out;
  };
  const at = (data, b) => X(shellMid(data, b));
  if (cover) {
    c.setLineDash([3, 3]);
    c.lineWidth = 1;
    c.globalAlpha = 0.55;
    for (const { data, color } of layers) {
      c.strokeStyle = color;
      c.beginPath();
      data.coverage.forEach((f, b) => {
        const y = C(f);
        if (b) c.lineTo(X(data.edges[b]), y);
        else c.moveTo(X(data.edges[b]), y);
        c.lineTo(X(data.edges[b + 1]), y);
      });
      c.stroke();
    }
    c.setLineDash([]);
    c.globalAlpha = 1;
  }
  for (const { data, color } of layers) {
    if (!band) break;
    c.fillStyle = color;
    c.globalAlpha = 0.16;
    for (const run of runs(data)) {
      if (!run.some((b) => data.sigma[b] > 0)) continue;
      c.beginPath();
      run.forEach((b, i) => {
        const y = Y(data.intensity[b] + (data.sigma[b] || 0));
        if (i) c.lineTo(at(data, b), y);
        else c.moveTo(at(data, b), y);
      });
      for (let i = run.length - 1; i >= 0; i--) {
        const b = run[i], v = data.intensity[b] - (data.sigma[b] || 0);
        c.lineTo(at(data, b), log && v <= 0 ? bottom + 10 : Y(v));
      }
      c.closePath();
      c.fill();
    }
    c.globalAlpha = 1;
  }
  // Markers when the bins are far enough apart to see them.
  const visible = layers[0].data.edges.filter((e, b, all) => b + 1 < all.length && (e + all[b + 1]) / 2 >= x0 && (e + all[b + 1]) / 2 <= x1).length;
  const markers = visible < aw / 6;
  for (const { data, color } of layers) {
    c.strokeStyle = c.fillStyle = color;
    c.lineWidth = 1.6;
    c.lineJoin = 'round';
    for (const run of runs(data)) {
      c.beginPath();
      run.forEach((b, i) => (i ? c.lineTo(at(data, b), Y(data.intensity[b])) : c.moveTo(at(data, b), Y(data.intensity[b]))));
      c.stroke();
      if (!markers && run.length > 1) continue;
      for (const b of run) {
        c.beginPath();
        c.arc(at(data, b), Y(data.intensity[b]), 2.2, 0, 2 * Math.PI);
        c.fill();
      }
    }
  }
  if (!exporting && plot.at !== null) {
    // The bin under the cursor.
    const x = at(layers[0].data, plot.at);
    c.strokeStyle = 'rgba(18, 24, 33, 0.35)';
    c.lineWidth = 1;
    c.setLineDash([4, 4]);
    c.beginPath(); c.moveTo(x, pad.t); c.lineTo(x, bottom); c.stroke();
    c.setLineDash([]);
    for (const { data, color } of layers) {
      const v = data.intensity[plot.at];
      if (!Number.isFinite(v) || (log && v <= 0)) continue;
      c.fillStyle = '#ffffff';
      c.strokeStyle = color;
      c.lineWidth = 2;
      c.beginPath(); c.arc(x, Y(v), 4, 0, 2 * Math.PI); c.fill(); c.stroke();
    }
  }
  const d = plot.drag;
  if (!exporting && d?.moved && clickMode !== 'move') {
    // The zoom box being dragged; a flat one zooms x only and spans the height.
    const flat = Math.abs(d.py1 - d.py) < 12;
    const xa = Math.min(d.px, d.px1), ya = flat ? pad.t : Math.min(d.py, d.py1);
    const bw = Math.abs(d.px1 - d.px), bh = flat ? ah : Math.abs(d.py1 - d.py);
    c.fillStyle = 'rgba(47, 116, 230, 0.12)';
    c.fillRect(xa, ya, bw, bh);
    c.setLineDash([5, 3]);
    c.strokeStyle = '#2f74e6';
    c.lineWidth = 1.5;
    c.strokeRect(xa, ya, bw, bh);
    c.setLineDash([]);
  }
  c.restore();

  // Axes, ticks and labels.
  c.strokeStyle = AXIS;
  c.lineWidth = 1;
  c.beginPath();
  c.moveTo(pad.l + 0.5, pad.t);
  c.lineTo(pad.l + 0.5, bottom + 0.5);
  c.lineTo(pad.l + aw, bottom + 0.5);
  if (cover) c.lineTo(pad.l + aw + 0.5, pad.t);
  c.stroke();
  c.fillStyle = INK2;
  c.font = `10.5px ${MONO}`;
  c.textAlign = 'center';
  c.textBaseline = 'top';
  c.beginPath();
  for (const t of xt) {
    const x = Math.round(X(t)) + 0.5;
    c.moveTo(x, bottom); c.lineTo(x, bottom + 4);
    c.fillText(fmt(t), x, bottom + 7);
  }
  c.textAlign = 'right';
  c.textBaseline = 'middle';
  for (const t of yt) {
    const y = Math.round(T(t)) + 0.5;
    c.moveTo(pad.l, y); c.lineTo(pad.l - 4, y);
    c.fillText(fmtValue(log ? 10 ** t : t), pad.l - 7, y);
  }
  if (cover) {
    c.textAlign = 'left';
    for (const f of [0, 0.5, 1]) {
      const y = Math.round(C(f)) + 0.5;
      c.moveTo(pad.l + aw, y); c.lineTo(pad.l + aw + 4, y);
      c.fillText(`${100 * f}%`, pad.l + aw + 7, y);
    }
  }
  c.stroke();
  c.fillStyle = INK;
  c.font = `600 11.5px ${SANS}`;
  c.textAlign = 'center';
  c.textBaseline = 'alphabetic';
  c.fillText(plot.xLabel(), pad.l + aw / 2, bottom + 36);
  const vertical = (text, x) => {
    c.save();
    c.translate(x, pad.t + ah / 2);
    c.rotate(-Math.PI / 2);
    c.fillText(text, 0, 0);
    c.restore();
  };
  vertical(plot.yLabel(log), 16);
  if (cover) {
    c.fillStyle = INK2;
    c.font = `11px ${SANS}`;
    c.save();
    c.translate(w - 8, pad.t + ah / 2);
    c.rotate(Math.PI / 2);
    c.fillText('shell coverage', 0, 0);
    c.restore();
  }

  // Legend: when comparing (on plots large enough to spare the room), and on exports.
  if ((layers.length > 1 && aw >= 280 && ah >= 160) || exporting) {
    c.font = `600 11.5px ${SANS}`;
    const items = layers.map((l) => [l.color, `${layers.length > 1 ? `${l.letter}  ` : ''}${shortName(l.name, 32)}`]);
    const lw = Math.max(...items.map(([, t]) => c.measureText(t).width)) + 38, lh = 18 * items.length + 10;
    const lx = pad.l + aw - lw - 8, ly = pad.t + 8;
    c.fillStyle = 'rgba(255, 255, 255, 0.9)';
    c.strokeStyle = 'rgba(18, 24, 33, 0.18)';
    c.beginPath(); c.roundRect(lx, ly, lw, lh, 6); c.fill(); c.stroke();
    c.textAlign = 'left';
    c.textBaseline = 'middle';
    items.forEach(([color, text], i) => {
      const y = ly + 14 + 18 * i;
      c.strokeStyle = color;
      c.lineWidth = 2;
      c.beginPath(); c.moveTo(lx + 10, y); c.lineTo(lx + 26, y); c.stroke();
      c.fillStyle = INK2;
      c.fillText(text, lx + 32, y);
    });
  }

  if (!exporting) return;
  c.textAlign = 'left';
  c.textBaseline = 'alphabetic';
  c.fillStyle = INK;
  c.font = `650 14px ${SANS}`;
  c.fillText(plot.title, 10, 24);
  const titleWidth = c.measureText(plot.title).width;
  c.font = `11px ${MONO}`;
  c.fillStyle = INK2;
  c.fillText(plot.extras(layers[0].data).join(' · '), 20 + titleWidth, 24);
}

/** The pointer in plot coordinates: x and t (intensity, or log10 of it). */
function curvePoint(plot, e) {
  const v = plot.view;
  if (!v) return null;
  const [px, py] = plotLocal(plot, e);
  return {
    px, py,
    x: v.x0 + ((px - v.pad.l) / v.aw) * (v.x1 - v.x0),
    t: v.y0 + ((v.pad.t + v.ah - py) / v.ah) * (v.y1 - v.y0),
    inside: px >= v.pad.l && px <= v.pad.l + v.aw && py >= v.pad.t && py <= v.pad.t + v.ah,
  };
}

/** Mark the bin under the pointer and read it in the hover chip; null clears it. */
function hoverCurve(plot, e) {
  const pt = e && curvePoint(plot, e), layers = shownLayers(plot), first = layers[0]?.data;
  const b = pt?.inside && first ? shellAt(first.edges, pt.x) : -1;
  const at = b >= 0 && b < first.intensity.length ? b : null;
  plot.hover.hidden = at === null;
  if (at !== null) plot.hover.textContent = plot.readout(at, layers);
  if (plot.at === at) return;
  plot.at = at;
  drawPlot(plot);
  plot.onHover?.();
}

// Dragging on a plot: a box zooms (a flat one x only), and in Move mode the plot pans.
function startCurveDrag(plot, e) {
  const pt = curvePoint(plot, e);
  if (e.button !== 0 || !pt?.inside) return;
  plot.drag = { px: pt.px, py: pt.py, px1: pt.px, py1: pt.py, x: pt.x, view: { ...plot.view }, moved: false };
  try { plot.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic or already released pointer */ }
  if (clickMode === 'move') plot.canvas.classList.add('panning');
}

function moveCurveDrag(plot, e) {
  const d = plot.drag;
  if (!d) return;
  [d.px1, d.py1] = plotLocal(plot, e);
  d.moved ||= Math.hypot(d.px1 - d.px, d.py1 - d.py) > (clickMode === 'move' ? 1 : 5);
  if (!d.moved) return;
  if (clickMode !== 'move') { drawPlot(plot); return; }
  const v = d.view, sx = ((d.px1 - d.px) / v.aw) * (v.x1 - v.x0), st = ((d.py1 - d.py) / v.ah) * (v.y1 - v.y0);
  setCurveZoom(plot, { x: [v.x0 - sx, v.x1 - sx], y: plot.zoom?.y ? [v.y0 + st, v.y1 + st] : null });
}

function endCurveDrag(plot) {
  const d = plot.drag;
  plot.drag = null;
  plot.canvas.classList.remove('panning');
  if (!d || clickMode === 'move') return;
  // The width of the bin where the drag started: zooms keep at least two such bins across.
  const v = d.view, edges = shownLayers(plot)[0]?.data.edges ?? [0, 0], s = Math.max(0, shellAt(edges, d.x)), dx = edges[s + 1] - edges[s];
  if (!d.moved) {
    // A click in Zoom mode zooms into x 2×, after a pause so a double-click (full range) does not also zoom.
    if (clickMode !== 'zoom') return;
    clearTimeout(plot.clickTimer);
    plot.clickTimer = setTimeout(() => {
      const half = Math.max((v.x1 - v.x0) / 4, dx), start = Math.max(plot.xMin, d.x - half);
      setCurveZoom(plot, { x: [start, start + 2 * half], y: null });
    }, 250);
    return;
  }
  const x = (px) => v.x0 + ((px - v.pad.l) / v.aw) * (v.x1 - v.x0), t = (py) => v.y0 + ((v.pad.t + v.ah - py) / v.ah) * (v.y1 - v.y0);
  const range = [x(Math.min(d.px, d.px1)), x(Math.max(d.px, d.px1))];
  // At least two bins across.
  if (range[1] - range[0] < 2 * dx) { drawPlot(plot); return; }
  setCurveZoom(plot, { x: range, y: Math.abs(d.py1 - d.py) < 12 ? null : [t(Math.max(d.py, d.py1)), t(Math.min(d.py, d.py1))] });
}

function setCurveZoom(plot, zoom) {
  plot.zoom = zoom;
  plot.zoomReset.hidden = !zoom;
  drawPlot(plot);
}

/** Hover, zoom and pan on a plot's canvas, and its redraw on resizes. */
function wirePlot(plot) {
  plot.canvas.onpointerdown = (e) => startCurveDrag(plot, e);
  plot.canvas.onpointermove = (e) => { hoverCurve(plot, e); moveCurveDrag(plot, e); };
  plot.canvas.onpointerup = () => endCurveDrag(plot);
  plot.canvas.onpointercancel = () => { plot.drag = null; plot.canvas.classList.remove('panning'); drawPlot(plot); };
  plot.canvas.onpointerleave = () => hoverCurve(plot, null);
  plot.canvas.ondblclick = () => { clearTimeout(plot.clickTimer); setCurveZoom(plot, null); };
  plot.zoomReset.onclick = () => setCurveZoom(plot, null);
  resized.observe(plot.canvas);
}

/** Save a plot as a PNG at 3×, with its title and legend. */
function savePlotPNG(plot, name) {
  if (!shownLayers(plot).length) return;
  const out = document.createElement('canvas'), s = uiScale();
  drawCurves(plot, out, Math.max(plot.canvas.clientWidth / s, 560), Math.max(plot.canvas.clientHeight / s, 360) + 32, 3, true);
  out.toBlob((blob) => download(blob, name));
}

// ---- Line cut -----------------------------------------------------------------------------------

// Directions a cut keeps to with Shift: small whole steps along the plane's axes.
const CUT_DIRECTIONS = [[1, 0], [0, 1], [1, 1], [1, -1], [1, 2], [2, 1], [1, -2], [2, -1]];
// A drag in Cut mode that starts this close (px) to an end of the cut moves that end.
const CUT_GRAB = 10;

/** The line cut view: its plot (see drawCurves()), header and footer. */
function setupCut() {
  const shell = viewShell('cut', slotSwitch('cut'), 'Line cut');
  shell.actions.innerHTML = `
    <button type="button" class="btn btn-ghost btn-xs zoom-reset" hidden title="Back to the full range (or double-click the plot)">Reset zoom</button>
    <div class="segmented sm cut-scale" role="group" aria-label="Intensity scale"><button type="button" data-value="linear" title="Linear intensity scale">Lin</button><button type="button" data-value="log" title="Logarithmic intensity scale">Log</button></div>
    <button type="button" class="icon-btn" title="Line cut options" aria-haspopup="dialog" aria-expanded="false" data-pop="pop-cut">${ICONS.gear}</button>
    <button type="button" class="icon-btn" title="Download the cut: CSV or text" aria-haspopup="menu" aria-expanded="false" data-pop="pop-cut-download">${ICONS.download}</button>
    <button type="button" class="icon-btn save" title="Save PNG">${ICONS.save}</button>`;
  shell.actions.append(shell.focusBtn, shell.maxBtn);
  shell.q('.view-body').innerHTML = `<canvas class="plot" role="img" aria-label="Line cut: the intensity along a line in a slice"></canvas><span class="overlay-chip" hidden></span>
    <button type="button" class="btn btn-primary view-cta" title="Switch the click mode to Cut, then drag across a slice">Draw a cut</button>`;
  const foot = document.createElement('footer');
  foot.className = 'view-foot';
  const labels = meta.dims.map((d) => d.label).join(', ');
  // Each label stays with its field when the footer wraps.
  foot.innerHTML = `
    <span class="foot-group"><span class="label">From</span>
      <input class="cut-from cut-end mono" type="text" spellcheck="false" autocomplete="off" aria-label="Cut start"
        title="Where the cut starts (${labels}). From and To share one coordinate: the slice the cut lies in"></span>
    <span class="foot-group"><span class="label">To</span>
      <input class="cut-to cut-end mono" type="text" spellcheck="false" autocomplete="off" aria-label="Cut end" title="Where the cut ends (${labels})"></span>
    <span class="foot-group grow"><span class="label" title="Diameter of the rod around the line whose voxels each point averages">W</span>
      <input class="cut-wslider foot-grow" type="range" min="0" max="1000" aria-label="Cut width (log scale)">
      <input class="cut-width num" type="number" min="0" step="any" aria-label="Cut width" title="Diameter of the rod around the line; empty for three voxels">
      <span class="unit cut-unit"></span></span>`;
  shell.section.append(foot);
  const q = shell.q;
  cut = {
    canvas: q('canvas'), hover: q('.overlay-chip'), caption: q('.view-caption'), zoomReset: q('.zoom-reset'), cta: q('.view-cta'),
    from: q('.cut-from'), to: q('.cut-to'), width: q('.cut-width'), wslider: q('.cut-wslider'), unit: q('.cut-unit'),
    a: newCurveLayer(), more: [newCurveLayer(), newCurveLayer(), newCurveLayer()], line: null, stale: false, download: false, zoom: null, drag: null, view: null, at: null,
    // As a plot (see drawCurves()):
    xMin: -Infinity, title: 'Line cut',
    log: () => cutScale === 'log', band: () => $('cut-band').checked, cover: () => false,
    end: (data) => data.edges[data.edges.length - 1],
    xLabel: () => { const { axis, path } = cutPath(); return `${withUnits(axis)} along ${path}`; },
    yLabel: (log) => (log ? 'Intensity, log scale' : 'Intensity'),
    empty: () => {
      if (!cut.line) return clickMode === 'cut' ? { text: 'Drag across a slice to draw a cut, or type its ends below.' } : null;
      const busy = [cut.a, ...cut.more].some((l) => l.busy), failed = layerOf(cut, shownIndex()).error;
      return { text: busy ? 'Computing the cut…' : failed || 'The cut appears here.', error: !busy && !!failed };
    },
    extras: (data) => {
      const { unit, scale } = cutWidthUnit();
      return [cutPath().path, `W ${fmt(cutSpec().radius * 2 * scale)} ${unit}`,
        data.order > 1 ? data.symmetry : 'no symmetry', data.masked && mask ? `mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : 'no mask'];
    },
    // The point: its coordinates, and I ± σ with the voxels pooled for each dataset.
    readout: (b, layers) => {
      const value = (d) => (Number.isFinite(d.intensity[b])
        ? `${fmtValue(d.intensity[b])}${Number.isFinite(d.sigma[b]) ? ` ± ${fmtValue(d.sigma[b])}` : ''} (${d.voxels[b]} vox)`
        : 'no data');
      const at = cutAt(shellMid(layers[0].data, b));
      return `(${at.map((x) => fmt(x)).join(', ')})  →  ${layers.map((l) => `${layers.length > 1 ? `${l.letter} ` : ''}${value(l.data)}`).join(' · ')}`;
    },
    // The point read here is marked on the slices the cut lies in.
    onHover: () => { for (const p of panels) if (cutLiesIn(p)) drawPanel(p); },
  };
  cut.cta.onclick = () => { setClickMode('cut'); persist(); };
  setSegmented(q('.cut-scale'), cutScale);
  segmented(q('.cut-scale'), (value) => {
    cutScale = value;
    if (cut.zoom) setCurveZoom(cut, cut.zoom.x ? { x: cut.zoom.x, y: null } : null);
    persist();
    drawPlot(cut);
  });
  // Typed ends and widths apply on Enter or leaving the field; the width slider as it moves.
  cut.from.onchange = cut.to.onchange = typeCutEnds;
  cut.width.onchange = () => { syncWidthSlider(); changeCutWidth(); };
  cut.wslider.oninput = () => {
    const [lo, hi] = widthRange();
    cut.width.value = sig(lo * (hi / lo) ** (Number(cut.wslider.value) / 1000), 2);
    changeCutWidth();
  };
  for (const input of [cut.from, cut.to, cut.width]) input.onkeydown = (e) => { if (e.key === 'Enter') input.blur(); };
  wirePlot(cut);
  q('.save').onclick = () => savePlotPNG(cut, `${shownStem()}_cut_${cutFileName()}.png`);
  showCut();
}

// cut.line = { a, b, key }: the cut's ends in display coordinates, a[dom] ≤ b[dom],
// and the slice it was drawn in (or typed ends lie in), or null for a cut free
// in 3-D. A cut in a slice lies in its plane and follows it.

/** The slice the cut was drawn in and follows, or null. */
const cutPanel = () => (cut?.line?.key ? panels.find((p) => p.key === cut.line.key) : null);

/** The cut's ends in display coordinates, those of a cut in a slice at the slice's current position. */
function cutEnds() {
  const { a, b } = cut.line, p = cutPanel();
  if (!p) return [a, b];
  const center = Number(p.center.value);
  return [a, b].map((end) => end.map((v, d) => (d === p.fixed ? center : v)));
}

/** Whether the cut lies in the plane of slice p. */
function cutLiesIn(p) {
  if (!cut?.line) return false;
  const center = Number(p.center.value);
  return cutEnds().every((end) => Math.abs(end[p.fixed] - center) < 1e-6);
}

/** A point (u, v) of slice p's plane in all three display coordinates. */
function inPlane(p, [u, v]) {
  const out = [0, 0, 0];
  out[p.fixed] = Number(p.center.value);
  out[p.x] = u;
  out[p.y] = v;
  return out;
}

/** Bin width of a display dimension (uniform bins). */
const binWidth = (dim) => (dim.edges[dim.edges.length - 1] - dim.edges[0]) / (dim.edges.length - 1);

/** A positive number typed in `input`, or null when it is empty (automatic). */
function typedPositive(input) {
  const text = input.value.trim(), x = Number(text);
  if (text === '') return null;
  if (!(x > 0)) throw new Error(`${input.getAttribute('aria-label')} must be a positive number, or empty for automatic.`);
  return x;
}

/**
 * How the cut's width W, the diameter of its rod, is given: in Å⁻¹ (with 2π)
 * with a lattice and in the axes' unit otherwise. `scale` turns lengths in
 * displayBasis() into that unit, and `auto` is three voxels across.
 */
function cutWidthUnit() {
  const { T, lattice } = displayBasis(), scale = lattice ? 2 * Math.PI : 1, units = meta.dims[0].units;
  const voxel = Math.max(...meta.dims.map((d, i) => binWidth(d) * Math.hypot(T[i], T[3 + i], T[6 + i])));
  return { T, scale, unit: lattice ? 'Å⁻¹' : meta.dims.every((d) => d.units === units) ? units || '' : '', auto: sig(3 * voxel * scale) };
}

/** The cut's geometry for lineCut() (see cut.js), in display coordinates; throws when its step or width cannot be read. */
function cutSpec() {
  const [a, b] = cutEnds(), dom = cutAxis(a, b);
  const { T, scale, auto } = cutWidthUnit();
  const step = typedPositive($('cut-step')) ?? binWidth(meta.dims[dom]);
  const width = typedPositive(cut.width) ?? auto;
  return { a, b, dom, edges: cutEdges(a, b, dom, step), radius: width / scale / 2, T };
}

/** The point of the cut's line, in display coordinates, where its axis has the value s. */
function cutAt(s) {
  const [a, b] = cutEnds(), dom = cutAxis(a, b), t = (s - a[dom]) / (b[dom] - a[dom]);
  return a.map((ad, d) => ad + t * (b[d] - ad));
}

const pointText = (pt) => pt.map((x) => String(Number(x.toFixed(4)))).join(', ');

/** The cut's axis, and its path as a function of that axis: (H, 0, 1), (H, H, H) or (H, 0.5H+1, 0). */
function cutPath() {
  const [a, b] = cutEnds(), dom = cutAxis(a, b), axis = meta.dims[dom];
  const term = (d) => {
    if (d === dom) return axis.label;
    const slope = (b[d] - a[d]) / (b[dom] - a[dom]), offset = a[d] - slope * a[dom];
    if (Math.abs(slope) < 1e-9) return fmt(a[d]);
    const k = Math.abs(slope - 1) < 1e-9 ? '' : Math.abs(slope + 1) < 1e-9 ? '−' : fmt(slope);
    return `${k}${axis.label}${Math.abs(offset) < 1e-9 ? '' : `${offset > 0 ? '+' : '−'}${fmt(Math.abs(offset))}`}`;
  };
  return { axis, path: `(${[0, 1, 2].map(term).join(', ')})` };
}

/** The cut's path for file names: H_0_1, H_-H+1_0. */
const cutFileName = () => cutPath().path.replace(/[()\s]/g, '').replace(/−/g, '-').replace(/,/g, '_');

/** Snap a point of slice p to the nearest voxel centre on each axis. */
function snapToVoxel(p, [u, v]) {
  const snap = (dim, x) => {
    const w = binWidth(dim), e0 = dim.edges[0];
    return Number((e0 + (Math.round((x - e0) / w - 0.5) + 0.5) * w).toFixed(6));
  };
  return [snap(meta.dims[p.x], u), snap(meta.dims[p.y], v)];
}

/**
 * The end of a cut from `from` toward `to` along the nearest of CUT_DIRECTIONS
 * (in the drawn geometry), a whole number of voxels along its axis.
 */
function snapDirection(p, from, to) {
  const g = geometry(p), sin = Math.sqrt(Math.max(0, 1 - g.cos * g.cos));
  const cart = ([u, v]) => [g.lx * u + g.ly * g.cos * v, g.ly * sin * v];
  const [dx, dy] = cart([to[0] - from[0], to[1] - from[1]]);
  let best = null;
  for (const dir of CUT_DIRECTIONS) {
    const [cx, cy] = cart(dir), t = (dx * cx + dy * cy) / (cx * cx + cy * cy), off = Math.hypot(dx - t * cx, dy - t * cy);
    if (!best || off < best.off) best = { dir, t, off };
  }
  const { dir, t } = best, dom = cutAxis([0, 0], dir), w = binWidth(meta.dims[dom ? p.y : p.x]) / Math.abs(dir[dom]);
  const k = Math.round(t / w) * w;
  return [from[0] + k * dir[0], from[1] + k * dir[1]].map((x) => Number(x.toFixed(6)));
}

/** Make the cut run between points P and Q of slice p's plane (either way round), in that slice. */
const setCutLine = (p, P, Q) => setCutEnds(inPlane(p, P), inPlane(p, Q), p);

/**
 * Make the cut run between display points P and Q (either way round), in
 * slice p or, without one, free in 3-D. The first cut shows in the fourth
 * place; later ones leave it as it is, so a cut can be changed while watching
 * the 3-D view.
 */
function setCutEnds(P, Q, p = null) {
  if (P.every((v, d) => v === Q[d])) return;
  const dom = cutAxis(P, Q), [a, b] = P[dom] <= Q[dom] ? [P, Q] : [Q, P];
  const first = !cut.line;
  cut.line = { key: p?.key ?? null, a, b };
  if (first && slot !== 'cut') setSlot('cut');
  syncWidthSlider();
  requestCut();
  // Every slice shows the cut, or where it crosses.
  panels.forEach(drawPanel);
  show3DCut();
}

// Cut mode on a slice: a drag from an end of the cut moves that end, and one
// from elsewhere draws a new cut. Ends snap to voxel centres, or with Shift to
// a lattice direction from the other end.
function startCutDrag(p, e) {
  if (e.button !== 0 || !p.inverse) return;
  const [x, y] = localPoint(p, e), ends = cutLiesIn(p) ? cutEnds().map((end) => [end[p.x], end[p.y]]) : [];
  const grabbed = ends.findIndex((end) => { const [ex, ey] = p.project(...end); return Math.hypot(ex - x, ey - y) < CUT_GRAB; });
  p.cutDrag = { x0: x, y0: y, moved: false, fresh: grabbed < 0, anchor: grabbed < 0 ? snapToVoxel(p, p.inverse(x, y)) : ends[1 - grabbed] };
  try { p.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic or already released pointer */ }
}

function moveCutDrag(p, e) {
  const d = p.cutDrag, [x, y] = localPoint(p, e);
  if (!d.moved) {
    if (Math.hypot(x - d.x0, y - d.y0) <= 5) return;
    d.moved = true;
    // A new cut starts at the full range.
    if (d.fresh) setCurveZoom(cut, null);
  }
  const to = p.inverse(x, y);
  setCutLine(p, d.anchor, e.shiftKey ? snapDirection(p, d.anchor, to) : snapToVoxel(p, to));
}

/**
 * Apply typed ends. When they share a coordinate, the cut lies in the slice
 * through it, which moves there; otherwise it is free in 3-D, like one from
 * (0, 0, 0) to (2, 2, 2) along (H, H, H).
 */
function typeCutEnds() {
  const read = (input) => {
    const values = input.value.replace(/[()[\]]/g, ' ').trim().split(/[\s,]+/).filter(Boolean).map(Number);
    if (values.length !== 3 || !values.every(Number.isFinite)) {
      throw new Error(`Write the cut's ends as three numbers (${meta.dims.map((d) => d.label).join(', ')}), like 1, 0, 0.`);
    }
    return values;
  };
  let P, Q;
  try {
    [P, Q] = [read(cut.from), read(cut.to)];
  } catch (err) {
    error(err.message);
    return;
  }
  error('');
  const p = [cutPanel(), ...panels].find((q) => q && Math.abs(P[q.fixed] - Q[q.fixed]) < 1e-9);
  if (p && Math.abs(Number(p.center.value) - P[p.fixed]) > 1e-9) {
    p.center.value = p.slider.value = P[p.fixed];
    paint(p.slider);
    request(p);
  }
  setCutEnds(P, Q, p ?? null);
}

// The width slider sets W on a log scale, from half a voxel to 30 voxels.
function widthRange() {
  const { auto } = cutWidthUnit();
  return [auto / 6, auto * 10];
}

function syncWidthSlider() {
  cut.wslider.disabled = !cut.line;
  if (!cut.line) return;
  const [lo, hi] = widthRange(), w = Number(cut.width.value) || cutWidthUnit().auto;
  cut.wslider.value = 1000 * clamp(Math.log(w / lo) / Math.log(hi / lo), 0, 1);
  paint(cut.wslider);
}

function changeCutWidth() {
  requestCut();
  panels.forEach(drawPanel);
  show3DCut();
}

/**
 * Recompute the cut after a change to its line, its slice, the processing or
 * its settings: now if its view is shown, otherwise once it is.
 */
function requestCut() {
  if (!cut) return;
  cut.stale = true;
  updateCut();
  showCut();
}

/** Send the cut requests when they are due: to A's worker, and the others' once their masks are built. */
function updateCut() {
  if (!cut?.stale || !cut.line || !plotShown(cut)) return;
  cut.stale = false;
  let spec;
  try {
    spec = cutSpec();
  } catch (err) {
    Object.assign(cut.a, { data: null, error: err.message });
    for (const layer of cut.more) layer.data = null;
    showCut();
    drawPlot(cut);
    return;
  }
  const removed = $('mask-removed').checked;
  queueCut(cut.a, worker, { cut: spec, maps: symmetry.maps, symmetry: symmetry.name, removed });
  for (const d of readyOthers()) {
    if (!d.maskBusy) queueCut(cut.more[d.slot], d.worker, { cut: spec, maps: d.maps, symmetry: d.maps.length > 1 ? symmetry.name : '1', removed });
  }
}

/** Ask a dataset's worker for the cut; one request in flight per dataset, the latest one next. */
function queueCut(layer, w, request) {
  layer.wanted = request;
  if (!layer.busy) sendCut(layer, w);
}

function sendCut(layer, w) {
  const q = layer.wanted;
  layer.wanted = null;
  layer.busy = true;
  w.postMessage({ type: 'cut', id: ++requestId, ...q });
}

/** The cut of dataset k (0 for A) has arrived, or failed with `message`. */
function cutResult(k, msg, message = '') {
  if (!cut) return;
  const layer = layerOf(cut, k);
  layer.busy = false;
  Object.assign(layer, msg ? { data: msg, error: '' } : { data: null, error: message });
  if (layer.wanted) sendCut(layer, workerOf(k));
  else endJob(`cut-${'abcd'[k]}`);
  showCut();
  drawPlot(cut);
  finishCutDownload();
}

/** The first cut of a file reads its uncertainties, which takes a moment on large files. */
function cutProgress(k, { label, fraction }) {
  if (cut && layerOf(cut, k).busy) setJob(`cut-${'abcd'[k]}`, k ? `Cut ${LETTERS[k]}` : 'Cut', label, fraction);
}

/** The cut view's header, fields and button. */
function showCut() {
  if (!cut) return;
  cut.cta.hidden = !!cut.line || clickMode === 'cut';
  syncWidthSlider();
  if (!cut.line) {
    cut.caption.textContent = clickMode === 'cut' ? 'drag across a slice' : '';
    cut.unit.textContent = '';
    return;
  }
  const [a, b] = cutEnds(), { unit, auto, scale } = cutWidthUnit();
  cut.unit.textContent = unit;
  cut.width.placeholder = String(auto);
  $('cut-step').placeholder = String(sig(binWidth(meta.dims[cutAxis(a, b)]), 3));
  if (document.activeElement !== cut.from) cut.from.value = pointText(a);
  if (document.activeElement !== cut.to) cut.to.value = pointText(b);
  const d = cut.a.data, { path } = cutPath();
  let width = auto;
  try { width = cutSpec().radius * 2 * scale; } catch { /* the field says why */ }
  if ([cut.a, ...cut.more].some((l) => l.busy)) {
    cut.caption.textContent = 'computing…';
  } else if (cut.a.error || !d) {
    cut.caption.textContent = cut.a.error || '';
  } else {
    cut.caption.textContent = `${path} · W ${fmt(width)} ${unit}`
      + `${d.order > 1 ? ` · ${d.symmetry}` : ''}${d.masked ? ' · masked' : ''}${d.removed ? ' · removed only' : ''}${d.errors ? ' · ±σ' : ''}`;
  }
  cut.caption.title = `Along ${path} from (${cut.from.value}) to (${cut.to.value}), averaging a rod ${fmt(width)} ${unit} across`
    + `${d ? `; ${d.intensity.length} points in ${d.seconds.toFixed(2)} s` : ''}`
    + readyOthers().filter((o) => cut.more[o.slot].error).map((o) => `\n${o.letter}: ${cut.more[o.slot].error}`).join('');
}

/**
 * The cut on slice p: in a slice whose plane it lies in, the band where its rod
 * meets the plane, its line and ends, and the point read on its plot; in the
 * others, where its rod passes through: the section and the point on its line.
 */
function drawCutOverlay(c, p, corners) {
  if (!cut?.line || !(plotShown(cut) || clickMode === 'cut')) return;
  let spec;
  try { spec = cutSpec(); } catch { return; }
  const flat = (pt) => p.project(pt[p.x], pt[p.y]);
  const outline = (points) => {
    c.beginPath();
    points.map(flat).forEach(([x, y], i) => (i ? c.lineTo(x, y) : c.moveTo(x, y)));
    c.closePath();
    c.fillStyle = 'rgba(255, 255, 255, 0.16)';
    c.fill();
    c.setLineDash([4, 3]);
    c.lineWidth = 1;
    c.strokeStyle = 'rgba(255, 255, 255, 0.85)';
    c.stroke();
    c.setLineDash([]);
  };
  const dot = ([x, y], r) => {
    c.beginPath(); c.arc(x, y, r, 0, 2 * Math.PI);
    c.fillStyle = '#ffffff'; c.fill();
    c.lineWidth = 1.5; c.strokeStyle = INK; c.stroke();
  };
  c.save();
  c.beginPath();
  corners.forEach(([x, y], i) => (i ? c.lineTo(x, y) : c.moveTo(x, y)));
  c.closePath();
  c.clip();
  if (!cutLiesIn(p)) {
    const crossing = cutCrossing(spec, p.fixed, Number(p.center.value));
    if (crossing) {
      outline(crossing.ring);
      dot(flat(crossing.at), 3.5);
    }
    c.restore();
    return;
  }
  const [A, B] = cutEnds().map(flat);
  outline(cutBand(spec, p.fixed));
  // A white line on a dark one, so it shows on dark and bright colors alike.
  c.lineCap = 'round';
  for (const [color, width] of [['rgba(18, 24, 33, 0.6)', 4], ['#ffffff', 2]]) {
    c.strokeStyle = color;
    c.lineWidth = width;
    c.beginPath(); c.moveTo(...A); c.lineTo(...B); c.stroke();
  }
  c.restore();
  c.save();
  for (const end of [A, B]) dot(end, 4.5);
  const first = cut.at !== null && shownLayers(cut)[0];
  if (first) {
    const [x, y] = flat(cutAt(shellMid(first.data, cut.at)));
    c.beginPath(); c.arc(x, y, 5.5, 0, 2 * Math.PI);
    c.lineWidth = 2.5; c.strokeStyle = '#ffffff'; c.stroke();
    c.lineWidth = 1.5; c.strokeStyle = first.color; c.stroke();
  }
  c.restore();
}

/** Save the cut as 'csv' or 'dat' (text), once it is computed. */
function downloadCut(format) {
  if (!cut?.line) return;
  cut.download = format;
  updateCut();
  finishCutDownload();
}

function finishCutDownload() {
  if (!cut?.download || cut.stale || [cut.a, ...readyOthers().map((d) => cut.more[d.slot])].some((l) => l.busy || l.wanted)) return;
  const format = cut.download;
  cut.download = false;
  if (!cut.a.data) return;
  const name = `${setStem(curveSets(cut).map(([, file]) => stemOf(file)))}_cut_${cutFileName()}`;
  if (format === 'csv') download(new Blob([cutCSV()], { type: 'text/csv' }), `${name}.csv`);
  else download(new Blob([cutText()], { type: 'text/plain' }), `${name}.dat`);
}

/**
 * The cut's points: per point its position along the cut's axis, its three
 * coordinates, and I, σ and voxels of each dataset shown (named in `values`).
 */
function cutTable() {
  const sets = curveSets(cut);
  const data = sets[0][2];
  // Positions are sums of steps: rounding leaves 1e-16 where 0 is meant.
  const clean = (x) => (Math.abs(x) < 1e-9 ? 0 : x);
  const rows = Array.from(data.intensity, (_, k) => {
    const s = shellMid(data, k);
    return [clean(s), ...cutAt(s).map(clean), ...sets.flatMap(([, , d]) => [d.intensity[k], d.sigma[k], d.voxels[k]])];
  });
  const values = sets.length > 1 ? sets.flatMap(([k]) => [`I_${k}`, `sigma_${k}`, `voxels_${k}`]) : ['I', 'sigma', 'voxels'];
  return { sets, data, rows, values };
}

/**
 * The cut as CSV: a header row, then per point its coordinates and I, σ and
 * voxels of each dataset, with empty cells where there is no value. The
 * position along the cut is the coordinate of its axis, so it is not repeated.
 */
function cutCSV() {
  const { rows, values } = cutTable();
  const cell = (x) => (Number.isFinite(x) ? String(Number(x.toPrecision(7))) : '');
  const quote = (text) => (/[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text);
  const lines = [[...meta.dims.map((d) => d.label), ...values].map(quote).join(',')];
  for (const row of rows) lines.push(row.slice(1).map(cell).join(','));
  return `${lines.join('\r\n')}\r\n`;
}

/** The cut as text: a commented header, then per point its position, coordinates, and I, σ and voxels of each dataset. */
function cutText() {
  const { sets, data, rows, values } = cutTable();
  const both = sets.length > 1, p = cutPanel(), spec = cutSpec(), { axis, path } = cutPath(), { unit, scale } = cutWidthUnit();
  const step = (data.edges[data.edges.length - 1] - data.edges[0]) / data.intensity.length;
  const where = p ? `drawn in the ${meta.dims[p.x].label}-${meta.dims[p.y].label} slice at ${meta.dims[p.fixed].label} = ${String(Number(spec.a[p.fixed].toPrecision(7)))}` : 'in 3-D';
  const num = (x) => (Number.isFinite(x) ? String(Number(x.toPrecision(7))) : 'nan');
  const lines = [
    `# Line cut along ${path} of ${listText(sets.map(([k, name]) => (both ? `${k} = ${name}` : name)))}`,
    `# Written by NeXus Viewer (https://drthyang.github.io/neutron-nexus-viewer/) on ${new Date().toISOString()}`,
    `# From (${pointText(spec.a)}) to (${pointText(spec.b)}), ${where}`,
    `# Points every ${num(step)} in ${axis.label}; each averages the voxels whose centres lie in a rod of diameter ${num(2 * spec.radius * scale)} ${unit || 'units'}`,
    `#   around the line (within ${num(spec.radius * scale)} of it, in any direction across it) and project within ±${num(step / 2)} ${axis.label} of the point.`,
    `# Symmetry: ${symmetry.ops.length > 1 ? `${symmetry.name} (${symmetry.ops.length} operations), equivalent voxels pooled` : 'none'}`,
    `# Mask: ${mask ? `coverage-edge erosion ${mask.radius}${mask.k ? `, outlier cut ${mask.k} sigma` : ''}${both ? ', built for each dataset' : ''}` : 'none'}${data.removed ? ' (removed voxels only)' : ''}`,
    '# I: the equal-weight mean of the distinct measured, unmasked voxels in the union of the symmetry orbits of the voxels at the point, as in the slices',
    `# sigma: ${sets.map(([k, , d]) => `${both ? `${k} ` : ''}${d.errors ? `propagated from ${d.errors}, voxels independent` : 'nan, no uncertainties in the file'}`).join('; ')}`,
    '# voxels: the number of distinct voxels with data pooled at the point',
    `# ${axis.label} ${meta.dims.map((d) => d.label).join(' ')} ${values.join(' ')}`,
  ];
  for (const row of rows) lines.push(row.map(num).join(' '));
  return `${lines.join('\n')}\n`;
}

// ---- Comparing datasets -------------------------------------------------------------------

const freeSlot = () => slots.indexOf(null);

/**
 * Open a file as dataset B, C or D (`slot` 0, 1 or 2, by default the first free
 * one; a dataset there is replaced) in its own worker; its slices follow A's
 * positions and processing.
 */
function openCompare(file, slot = freeSlot()) {
  if (!panels.length || slot < 0) return; // viewer not ready yet, or no slot free
  closeCompare(slot);
  const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  const d = slots[slot] = {
    slot, letter: LETTERS[slot + 1], worker: w, name: file.name, size: file.size, meta: null, ready: false, autoscaled: false,
    maps: [IDENTITY_MAP], mapsNote: '', mask: null, maskNote: '', maskBusy: false, work: null,
  };
  const on = compareHandlers(d);
  w.onmessage = ({ data }) => { if (slots[slot] === d) on[data.type]?.(data); };
  w.onerror = (e) => {
    e.preventDefault();
    if (slots[slot] !== d) return;
    closeCompare(slot);
    error(`Dataset ${d.letter}: ${e.message || 'the HDF5 reader could not start.'}`);
  };
  // Axes that are A's in another order are displayed in A's order.
  w.postMessage({ type: 'open', file, order: meta.dims.map((x) => x.label) });
  compareWork(d, 'open', 'Opening', 0);
}

/** Work on dataset d (`kind` 'open' or 'mask'): a job in the top bar, and the step under its name. */
function compareWork(d, kind, step, fraction) {
  setJob(`${kind}-${d.letter.toLowerCase()}`, `${kind === 'mask' ? 'Mask' : 'Opening'} ${d.letter}`, step, fraction);
  d.work = { step, fraction };
  showCompare();
}

async function openCompareURL(url, slot = freeSlot()) {
  if (!panels.length || slot < 0) return;
  closeCompare(slot);
  const d = slots[slot] = { slot, letter: LETTERS[slot + 1], worker: null, name: urlName(url), ready: false, work: null };
  compareWork(d, 'open', 'Downloading', 0);
  try {
    const file = await fetchFile(url, (got, total) => { if (slots[slot] === d) compareWork(d, 'open', `Downloading ${downloaded(got, total)}`, total ? got / total : 0); });
    if (slots[slot] === d) openCompare(file, slot);
  } catch (err) {
    if (slots[slot] !== d) return;
    closeCompare(slot);
    error(`Could not download ${url}: ${err.message}. The server must allow cross-origin requests.`);
  }
}

/** Remove the dataset in `slot` (B, C or D), if any. */
function closeCompare(slot) {
  slots[slot]?.worker?.terminate();
  slots[slot] = null;
  const l = 'bcd'[slot];
  if (compareView === l) compareView = 'split';
  endJob(`open-${l}`, `mask-${l}`, `iq-${l}`, `cut-${l}`);
  for (const p of panels) {
    p.more[slot] = newLayer();
    showCaption(p);
  }
  showCompare();
  if (meta && panels.length) {
    powder.more[slot] = newCurveLayer();
    cut.more[slot] = newCurveLayer();
    showCut();
    requestPowder();
    describe();
    redraw();
  }
}

/** The symmetry's index maps on dataset d's grid, or the identity (with a note) when they do not fit it. */
function compareMaps(d) {
  if (symmetry.ops.length === 1) return { maps: [IDENTITY_MAP], mapsNote: '' };
  try {
    return { maps: indexMaps(symmetry.ops, d.meta.dims), mapsNote: '' };
  } catch (err) {
    return { maps: [IDENTITY_MAP], mapsNote: `Symmetry not applied to ${d.letter}: ${err.message}` };
  }
}

function sendCompareMask(d, radius, k) {
  const kd = d.maps.length >= 3 ? k : 0;
  d.maskNote = k && !kd ? `Outlier cut not applied to ${d.letter}: it needs a symmetry of at least 3 operations that fits ${d.letter}'s grid.` : '';
  d.maskBusy = true;
  d.maskRequested = `${radius},${k}`;
  d.worker.postMessage({ type: 'mask', id: ++requestId, radius, k: kd, maps: d.maps, symmetry: symmetry.name });
  compareWork(d, 'mask', 'Building mask', 0);
}

/** The handlers of dataset d's worker messages. */
function compareHandlers(d) {
  const k = d.slot + 1, l = d.letter.toLowerCase();
  const done = (kind) => {
    d.work = null;
    endJob(`${kind}-${l}`);
  };
  return {
    progress: ({ label, fraction }) => compareWork(d, 'open', label, fraction),
    meta: ({ info }) => { d.meta = info; },
    ready: ({ stats, seconds }) => {
      Object.assign(d.meta, { stats, seconds });
      d.ready = true;
      done('open');
      Object.assign(d, compareMaps(d));
      for (const p of panels) p.more[d.slot] = newLayer();
      // It gets A's mask before its first slices.
      if (mask) sendCompareMask(d, mask.radius, mask.k);
      else panels.forEach((p) => request(p, d));
      showCompare();
      powder.more[d.slot] = newCurveLayer();
      cut.more[d.slot] = newCurveLayer();
      requestPowder();
      requestCut();
      describe();
      redraw();
    },
    slice: (msg) => {
      const p = panels.find((q) => q.fixed === msg.fixed);
      if (!p) return;
      const layer = p.more[d.slot];
      Object.assign(layer, { busy: false, data: msg, error: '', version: layer.version + 1 });
      showCaption(p);
      if (layer.wanted) sendOther(p, d);
      // Once its first slices are in, the automatic range covers it too.
      if (!d.autoscaled && panels.every((q) => q.more[d.slot].data || q.more[d.slot].error)) {
        d.autoscaled = true;
        if (rangeIsAuto) autoRange();
      }
      redraw();
    },
    error: ({ fixed, message }) => {
      const p = panels.find((q) => q.fixed === fixed);
      if (!p) {
        closeCompare(d.slot);
        error(`Dataset ${d.letter}: ${message}`);
        return;
      }
      const layer = p.more[d.slot];
      Object.assign(layer, { busy: false, data: null, error: message, version: layer.version + 1 });
      showCaption(p);
      if (layer.wanted) sendOther(p, d);
      redraw();
    },
    'progress-mask': ({ label, fraction }) => compareWork(d, 'mask', label, fraction),
    mask: ({ stats, radius, k: cutK }) => {
      d.maskBusy = false;
      done('mask');
      d.mask = stats ? { ...stats, radius, k: cutK } : null;
      d.maskVersion = (d.maskVersion ?? 0) + 1;
      showCompare();
      describe();
      panels.forEach((p) => request(p, d));
      requestPowder();
      requestCut();
    },
    'mask-error': ({ message }) => {
      d.maskBusy = false;
      done('mask');
      d.maskNote = `Mask for ${d.letter} failed: ${message}`;
      showCompare();
      panels.forEach((p) => request(p, d));
      requestPowder();
      requestCut();
    },
    'progress-powder': (msg) => powderProgress(k, msg),
    powder: (msg) => powderResult(k, msg),
    'powder-error': ({ message }) => powderResult(k, null, message),
    'progress-cut': (msg) => cutProgress(k, msg),
    cut: (msg) => cutResult(k, msg),
    'cut-error': ({ message }) => cutResult(k, null, message),
  };
}

// The chips of B, C and D in the top bar, made from the template in index.html.
const chips = [0, 1, 2].map((slot) => {
  const root = $('ds-template').content.firstElementChild.cloneNode(true);
  root.dataset.slot = slot;
  root.querySelector('.ab-tag').textContent = LETTERS[slot + 1];
  $('compare-open').before(root);
  const q = (sel) => root.querySelector(sel);
  return { root, main: q('.ds-main'), name: q('.ds-name'), facts: q('.ds-facts'), warn: q('.ds-warn'), close: q('.icon-btn') };
});

/**
 * The datasets compared with A in the top bar (each chip shows its name and
 * facts, or the work in progress on it), the Compare button that adds one, and
 * the A / Split / B … control.
 */
function showCompare() {
  const list = others(), ready = readyOthers(), quadrants = list.length > 1;
  document.body.classList.toggle('comparing', ready.length > 0);
  document.body.classList.toggle('comparing-more', list.length > 1);
  $('compare-open').hidden = list.length === slots.length || !panels.length;
  $('compare-open').title = `Compare with ${list.length ? 'another' : 'a second'} file (${LETTERS[freeSlot() + 1] ?? 'D'}): two datasets split every slice along its diagonal, `
    + 'three or four share it in quadrants, and I(Q) and the cut show them all. Or drop a file here.';
  $('compare-view').hidden = $('file-sep').hidden = !ready.length;
  $('file-sep').innerHTML = quadrants ? ICONS.quadrants : ICONS.diagonal;
  $('file-sep').title = quadrants ? 'Compared: each slice is shared in quadrants, A, B, C, D from lower left to upper right' : 'Compared: each slice is split along its diagonal';
  for (const button of $('compare-view').querySelectorAll('button')) {
    const k = 'abcd'.indexOf(button.dataset.value);
    if (k > 0) button.hidden = !slots[k - 1]?.ready;
  }
  const splitButton = $('compare-view').querySelector('[data-value="split"]');
  splitButton.title = quadrants
    ? 'Share each slice in quadrants: A lower left, B lower right, C upper left, D upper right'
    : 'Split each slice along its diagonal: A below, the other dataset above';
  splitButton.querySelector('.split-icon').innerHTML = quadrants ? ICONS.quadrants : ICONS.diagonal;
  setSegmented($('compare-view'), compareShown() ?? 'split');
  $('dataset-tag').hidden = !list.length;
  if (sourceName) setFileName($('dataset-name'), sourceName, list[0]?.name ?? null);
  const title3d = views['3d']?.section.querySelector('.view-title');
  if (title3d) title3d.textContent = ready.length ? 'Isosurface · A' : 'Isosurface';
  const labels = (dims) => dims.map((x) => x.label).join(' ');
  slots.forEach((d, slot) => {
    const chip = chips[slot];
    chip.root.hidden = !d;
    if (!d) return;
    setFileName(chip.name, d.name, sourceName);
    const step = d.work && `${d.work.step}${d.work.fraction > 0 ? ` · ${Math.round(100 * d.work.fraction)}%` : '…'}`;
    chip.facts.textContent = step || (d.ready ? datasetFacts(d.meta) : 'Opening…');
    chip.close.title = d.ready ? `Remove ${d.letter}` : 'Cancel';
    // Warnings: its axes, and processing that could not be applied to it.
    const warnings = [];
    if (d.ready && labels(d.meta.dims) !== labels(meta.dims)) warnings.push(`${d.letter}'s axes (${labels(d.meta.dims)}) differ from A's (${labels(meta.dims)}); ${d.letter} is drawn on A's axes.`);
    if (d.mapsNote) warnings.push(d.mapsNote);
    if (d.maskNote) warnings.push(d.maskNote);
    chip.warn.hidden = !warnings.length;
    chip.warn.title = warnings.join('\n');
    chip.main.title = `Dataset ${d.letter}: ${d.name}${d.size ? ` (${mb(d.size)})` : ''}`
      + `${d.ready ? `\n${datasetDetails(d.meta, d.mask)}` : ''}${warnings.length ? `\n${warnings.join('\n')}` : ''}`
      + `\nSlice positions, symmetry, mask settings and the color scale apply to all datasets.\nClick to replace ${d.letter} with another file.`;
  });
}

// ---- Startup -------------------------------------------------------------------------

function restore() {
  for (const name of Object.keys(LUTS)) $('cmap').append(new Option(name, name));
  for (const [name] of PRESETS) $('sym-preset').append(new Option(name === '1' ? 'None (1)' : name, name));
  $('sym-preset').append(new Option('Custom', 'custom'));
  let saved = {};
  try { saved = JSON.parse(localStorage.getItem('nxv-settings')) ?? {}; } catch { /* storage unavailable */ }
  if (['asinh', 'linear', 'log'].includes(saved.scale)) setSegmented($('scale'), saved.scale);
  if (saved.cmap in LUTS) $('cmap').value = saved.cmap;
  cmapChoice = $('cmap').value;
  if (typeof saved.angles === 'boolean') $('angles').checked = saved.angles;
  if (typeof saved.guides === 'boolean') $('guides').checked = saved.guides;
  if (typeof saved.grid === 'boolean') $('grid').checked = saved.grid;
  setClickMode(['zoom', 'move', 'cut'].includes(saved.clickMode) ? saved.clickMode : 'navigate');
  if (['linear', 'log'].includes(saved.iqScale)) powderScale = saved.iqScale;
  if (['1', '2', '3'].includes(saved.iqSplit)) setSegmented($('iq-split'), saved.iqSplit);
  if (typeof saved.iqBand === 'boolean') $('iq-band').checked = saved.iqBand;
  if (typeof saved.iqCoverage === 'boolean') $('iq-coverage').checked = saved.iqCoverage;
  if (['linear', 'log'].includes(saved.cutScale)) cutScale = saved.cutScale;
  if (typeof saved.cutBand === 'boolean') $('cut-band').checked = saved.cutBand;
  try {
    const l = JSON.parse(localStorage.getItem('nxv-layout'));
    if (['quad', 'focus'].includes(l?.mode)) lastMulti = layout = l.mode;
    if (SLOT_VIEWS.includes(l?.slot)) slot = l.slot;
    if (['hk', 'hl', 'kl', ...SLOT_VIEWS].includes(l?.key)) primary = l.key;
  } catch { /* storage unavailable: keep the defaults */ }
  paintColorbar();
}

function persist() {
  const values = {
    cmap: cmapChoice, scale: $('scale').dataset.value, angles: $('angles').checked, guides: $('guides').checked, grid: $('grid').checked, clickMode,
    iqScale: powderScale, iqSplit: $('iq-split').dataset.value, iqBand: $('iq-band').checked, iqCoverage: $('iq-coverage').checked,
    cutScale, cutBand: $('cut-band').checked,
  };
  try { localStorage.setItem('nxv-settings', JSON.stringify(values)); } catch { /* storage unavailable */ }
}

restore();
show('intro');
status('', 'no file');
$('open').onclick = $('open-intro').onclick = () => $('file').click();
$('open-example').onclick = openDemo;
$('file').onchange = () => {
  if ($('file').files[0]) {
    pendingCompare = [];
    pendingProcessing = null;
    openFile($('file').files[0]);
  }
  $('file').value = '';
};
// The file chosen next opens as the dataset in `fileSlot`: added by Compare…, or replacing a chip's.
let fileSlot = 0;
const chooseFile = (slot) => {
  fileSlot = slot;
  $('file-b').click();
};
$('compare-open').onclick = () => chooseFile(freeSlot());
chips.forEach((chip, slot) => {
  chip.main.onclick = () => chooseFile(slot);
  chip.close.onclick = () => closeCompare(slot);
});
$('file-b').onchange = () => { if ($('file-b').files[0]) openCompare($('file-b').files[0], fileSlot); $('file-b').value = ''; };
segmented($('compare-view'), (value) => {
  compareView = value;
  redraw();
});
$('error-close').onclick = () => error('');
for (const id of ['cmap', 'vmin', 'vmax', 'soft', 'limit', 'angles', 'guides', 'grid']) $(id).addEventListener('input', redraw);
$('cmap').addEventListener('change', () => { cmapChoice = $('cmap').value; });
for (const id of ['cmap', 'angles', 'guides', 'grid']) $(id).addEventListener('change', persist);
for (const id of ['vmin', 'vmax', 'soft']) $(id).addEventListener('input', () => { rangeIsAuto = false; });
$('angles').addEventListener('change', () => { if (meta) { describe(); requestCut(); } });
$('cmap').addEventListener('input', paintColorbar);
segmented($('click-mode'), (mode) => { setClickMode(mode); persist(); });
segmented($('layout'), (mode) => setLayout(mode, primary));
segmented($('scale'), (value) => {
  if (value === 'log' && !(Number($('vmin').value) > 0)) $('vmin').value = sig(Number($('vmax').value) / 1000 || 1);
  persist();
  redraw();
});
segmented($('iso-grid'), () => { if (iso) { iso.userLevel = false; requestIso(); } });
$('iso-slices').onchange = update3D;
$('auto').onclick = autoRange;
function setPanel(open) {
  document.body.classList.toggle('panel-collapsed', !open);
  $('panel-toggle').setAttribute('aria-pressed', String(open));
  $('panel-toggle').title = open ? 'Hide the control panel' : 'Show the control panel';
  try { localStorage.setItem('nxv-panel', open ? 'open' : 'closed'); } catch { /* storage unavailable */ }
}
$('panel-toggle').onclick = () => setPanel(document.body.classList.contains('panel-collapsed'));
// Collapsible panel sections, remembered per browser.
let collapsedSections = [];
try { collapsedSections = JSON.parse(localStorage.getItem('nxv-sections')) ?? []; } catch { /* storage unavailable */ }
for (const sec of document.querySelectorAll('.psec')) {
  const head = sec.querySelector('.psec-head');
  const set = (collapsed) => {
    sec.classList.toggle('collapsed', collapsed);
    head.setAttribute('aria-expanded', String(!collapsed));
  };
  set(collapsedSections.includes(sec.dataset.sec));
  head.onclick = () => {
    set(!sec.classList.contains('collapsed'));
    const now = [...document.querySelectorAll('.psec.collapsed')].map((s) => s.dataset.sec);
    try { localStorage.setItem('nxv-sections', JSON.stringify(now)); } catch { /* storage unavailable */ }
  };
}
try { if (localStorage.getItem('nxv-panel') === 'closed') setPanel(false); } catch { /* storage unavailable */ }
/** Show a panel element: open the panel and its section, scroll to it and flash it. */
function reveal(target) {
  setPanel(true);
  const sec = target.closest('.psec');
  if (sec?.classList.contains('collapsed')) sec.querySelector('.psec-head').click();
  target.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  target.classList.remove('flash');
  void target.offsetWidth;
  target.classList.add('flash');
}
for (const li of document.querySelectorAll('.pipeline li[data-target]')) li.onclick = () => reveal($(li.dataset.target));
$('sym-preset').onchange = () => {
  const preset = PRESETS.find(([name]) => name === $('sym-preset').value);
  if (!preset) { $('sym-ops').focus(); return; }
  $('sym-ops').value = preset[0] === '1' ? '' : preset[1];
  applySymmetry();
};
$('sym-ops').onkeydown = (e) => { if (e.key === 'Enter') applySymmetry(); };
$('sym-apply').onclick = applySymmetry;
$('mask-apply').onclick = () => applyMask();
$('mask-clear').onclick = () => applyMask(true);
for (const id of ['mask-erode', 'mask-k']) $(id).onkeydown = (e) => { if (e.key === 'Enter') applyMask(); };
$('mask-removed').onchange = () => { updateStates(); panels.forEach(request); };
$('mask-download').onclick = () => worker.postMessage({ type: 'mask-download' });
$('export-run').onclick = () => runExport();
$('export-open').onclick = openInNebula;
$('powder-show').onclick = showPowderView;
$('powder-download').onclick = downloadPowder;
segmented($('iq-split'), () => { persist(); requestPowder(); });
for (const id of ['iq-band', 'iq-coverage']) $(id).onchange = () => { persist(); drawPlot(powder); };
$('cut-band').onchange = () => { persist(); drawPlot(cut); };
for (const item of $('pop-cut-download').querySelectorAll('[data-format]')) {
  item.onclick = () => { closePopovers(); downloadCut(item.dataset.format); };
}
$('cut-step').onchange = changeCutWidth;
$('cut-step').onkeydown = (e) => { if (e.key === 'Enter') e.target.blur(); };

// Dropping a file opens it, or opens it as dataset B, C or D over that dataset's
// chip (replacing it) or over the Compare button (as the next free one).
const dropSlot = (e) => {
  if (!panels.length) return -1;
  const chip = e.target.closest?.('.ds-b');
  if (chip) return Number(chip.dataset.slot);
  return e.target.closest?.('#compare-open') ? freeSlot() : -1;
};
document.addEventListener('dragover', (e) => {
  e.preventDefault();
  const slot = dropSlot(e);
  document.body.classList.add('dragging');
  document.body.classList.toggle('drop-b', slot >= 0);
  if (slot >= 0) document.body.dataset.drop = LETTERS[slot + 1];
});
document.addEventListener('dragleave', (e) => { if (!e.relatedTarget) document.body.classList.remove('dragging', 'drop-b'); });
document.addEventListener('drop', (e) => {
  e.preventDefault();
  document.body.classList.remove('dragging', 'drop-b');
  const file = e.dataTransfer?.files?.[0], slot = dropSlot(e);
  if (!file) return;
  if (slot >= 0) openCompare(file, slot);
  else {
    pendingCompare = [];
    pendingProcessing = null;
    openFile(file);
  }
});

// ?url= opens a remote file, ?compare= (up to three times) the files to compare
// it with as datasets B, C and D, and ?sym= (a Laue class, such as 6/mmm) and
// ?mask= (erosion radius, optionally ",k" for the outlier cut) process them.
// ?demo opens the example.
const params = new URLSearchParams(location.search), remote = params.get('url');
if (params.has('demo')) openDemo();
else if (remote) {
  pendingCompare = params.getAll('compare').slice(0, slots.length);
  pendingProcessing = { sym: params.get('sym'), mask: params.get('mask') };
  openURL(remote);
}
