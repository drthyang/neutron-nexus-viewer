import { COLORMAPS } from './colormaps.js';
import { exportPlan } from './export.js';
import { cartesianBasis, nominalCell, planeGeometry, reciprocalMetric } from './nexus.js';
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
};

let worker = null, meta = null, panels = [], settings = null, sourceName = '', sourceSize = 0, autoscaled = false;
let requestId = 0, symmetry = NO_SYMMETRY, mask = null;
// Workspace layout: 'quad', 'focus' or 'single', around the primary view ('hk', 'hl', 'kl' or '3d').
let layout = 'quad', primary = 'hk', lastMulti = 'quad';
const views = {};
// What clicking a slice does: 'navigate' moves the other two slices, 'zoom' zooms in,
// 'move' drags the visible region.
let clickMode = 'navigate';
// 3-D view state: the lazily loaded View3D, its elements and the isosurface request queue.
let view3d = null, iso = null;
// Second dataset (B) for comparison, with its own worker, index maps and mask.
// `compareView` is what the slices show: 'split' (A below the diagonal, B above), 'a' or 'b'.
let compare = null, compareView = 'split', pendingCompare = null;
// Symmetry and mask to apply once the next file opens (?sym= and ?mask=, or the demo).
let pendingProcessing = null;
const DEMO = { url: 'examples/demo_300K.nxs', compare: 'examples/demo_10K.nxs', sym: '6/mmm', mask: '1' };
// True while the color range is the automatic one (not edited by hand).
let rangeIsAuto = false;
// The export for NEBULA3D in progress: { name }.
let exportJob = null;
const newLayer = () => ({ data: null, version: 0, busy: false, wanted: null, image: null, imageKey: null, error: '' });
const compareShown = () => (compare?.ready ? compareView : null);

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

function status(state, text) {
  $('status-dot').className = `status-dot ${state}`;
  $('status-text').textContent = text;
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
  status(meta ? 'ok' : '', meta ? 'ready' : 'no file');
  if (!panels.length) show('intro');
}

function openFile(file) {
  worker?.terminate();
  compare?.worker?.terminate();
  compare = null;
  exportJob = null;
  showCompare();
  closePopovers();
  meta = null;
  panels = [];
  autoscaled = false;
  symmetry = NO_SYMMETRY;
  mask = null;
  iso = null;
  sourceName = file.name;
  sourceSize = file.size;
  $('workspace').replaceChildren();
  error('');
  show('loading');
  status('busy', 'reading');
  setFileName($('dataset-name'), file.name);
  $('loading-title').textContent = 'Opening file';
  setFileName($('loading-name'), file.name);
  $('loading-size').textContent = mb(file.size);
  progress('Starting', 0);
  document.title = `${file.name} · NeXus Slice Viewer`;
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
    redraw();
    panels.forEach(request);
    setup3D();
    if (pendingCompare) openCompareURL(pendingCompare);
    pendingCompare = null;
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
    $('mask-status').className = 'note';
    $('mask-status').textContent = `${label}… ${Math.round(100 * fraction)}%`;
  },
  mask: ({ stats, radius, k, symmetry: group, seconds }) => {
    status('ok', 'ready');
    $('mask-apply').disabled = false;
    mask = stats ? { ...stats, radius, k, group } : null;
    // B may have opened while this mask was being built.
    if (compare?.ready && compare.maskRequested !== `${radius},${k}`) sendCompareMask(radius, k);
    $('mask-clear').disabled = $('mask-download').disabled = $('mask-removed').disabled = !mask;
    if (!mask) $('mask-removed').checked = false;
    showMask(seconds);
    panels.forEach((p) => request(p, 'a'));
    requestIso();
    describe();
  },
  'mask-error': ({ message }) => {
    status('ok', 'ready');
    $('mask-apply').disabled = false;
    $('mask-status').className = 'note error';
    $('mask-status').textContent = message;
    updateStates();
  },
  'mask-file': ({ blob }) => download(blob, `${stem()}_mask.npy.gz`),
  'progress-export': ({ label, fraction }) => showExportProgress(label, fraction),
  'export-file': ({ blob, stats, seconds }) => {
    const { name } = exportJob;
    exportJob = null;
    download(blob, name);
    status('ok', 'ready');
    showExport();
    $('export-state').textContent = 'saved';
    $('export-state').className = 'state ok';
    $('export-status').className = 'note ok';
    $('export-status').textContent = `Saved ${name} (${mb(blob.size)}, ${pct(stats.valid / stats.total)} of voxels valid) in ${seconds.toFixed(1)} s. `
      + 'In NEBULA3D, open it with Load volume…';
  },
  'export-error': ({ message }) => {
    exportJob = null;
    status('ok', 'ready');
    showExport();
    $('export-status').className = 'note error';
    $('export-status').textContent = message;
  },
};

/** Ask for the panel's slice of dataset 'a', 'b' or 'both' (B only when it is loaded). */
function request(p, which = 'both') {
  const center = Number(p.center.value), thickness = Number(p.width.value);
  if (!Number.isFinite(center) || !Number.isFinite(thickness) || thickness <= 0) {
    p.caption.textContent = 'Enter a finite center and a positive thickness.';
    return;
  }
  if (which !== 'b') {
    p.wanted = { center, thickness };
    if (!p.busy) send(p);
  }
  if (which !== 'a' && compare?.ready) {
    p.b.wanted = { center, thickness };
    if (!p.b.busy) sendB(p);
  }
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

function sendB(p) {
  const b = p.b, q = b.wanted;
  b.wanted = null;
  b.busy = true;
  compare.worker.postMessage({
    type: 'slice', id: ++requestId, fixed: p.fixed, ...q,
    maps: compare.maps, symmetry: compare.maps.length > 1 ? symmetry.name : '1', removed: $('mask-removed').checked,
  });
}

/** View header: position, bins and coverage, and B's coverage when comparing. */
function showCaption(p) {
  const F = meta.dims[p.fixed], d = p.data;
  if (!d) return;
  p.pos.textContent = `${F.label} = ${fmt(d.center)}`;
  const slab = (x) => `${F.label} ∈ [${fmt(x.slab[0])}, ${fmt(x.slab[1])}], ${x.bins} bin${x.bins === 1 ? '' : 's'}, ${pct(x.coverage)} of the plane measured`
    + `${x.order > 1 ? `, averaged over ${x.symmetry} (${x.order} operations)` : ''}${x.removed ? ', removed voxels only' : ''}`;
  let text = `${d.bins} bin${d.bins === 1 ? '' : 's'} · ${pct(d.coverage)}${d.order > 1 ? ` · ${d.symmetry}` : ''}${d.removed ? ' · removed only' : ''}`;
  let title = slab(d);
  if (compare?.ready) {
    const b = p.b.data;
    text = `A ${text} | B ${b ? pct(b.coverage) : p.b.error ? 'no data' : '…'}`;
    title = `A: ${title}\nB: ${b ? slab(b) : p.b.error || 'loading'}`;
  }
  p.caption.textContent = text;
  p.caption.title = title;
}

const symmetryNote = (data) => (data.order > 1 ? ` · ${data.symmetry} (${data.order} ops)` : '');
const stemOf = (name) => name.replace(/\.[^.]+$/, '');
const stem = () => stemOf(sourceName);

// ---- Viewer setup -------------------------------------------------------------

function describe() {
  const { dims, shape, lattice, stats } = meta;
  const bins = [0, 1, 2].map((d) => shape[2 - d]);
  const widths = dims.map((d) => fmt((d.edges[d.edges.length - 1] - d.edges[0]) / (d.edges.length - 1), 4));

  // Dataset section of the control panel: one table, and when comparing, a value
  // shared by A and B appears once while differing values get a line each.
  const { a, b, c, alpha, beta, gamma } = lattice ?? {};
  const recip = recipOf(lattice);
  const A = datasetFacts(meta, mask), B = compare?.ready ? datasetFacts(compare.meta, compare.mask) : null;
  const rows = [['Cell', 'cell'], ['Recip.', 'recip'], ['Grid', 'grid'], ['Measured', 'measured']];
  if (B && (mask || compare.mask)) rows.push(['Masked', 'masked']);
  $('data-kv').innerHTML = rows.map(([label, key]) => {
    const va = A[key] ?? '—', vb = B ? B[key] ?? '—' : null;
    let dd;
    if (!B) dd = `<dd>${va}${key === 'measured' ? ' of voxels' : ''}</dd>`;
    else if (va === vb) dd = `<dd title="Same for A and B">${va}</dd>`;
    else if (!/<br>/.test(va + vb) && (va + vb).length < 28) dd = `<dd><span class="ab-inline"><span class="ab-mini">A</span>${va}</span><span class="ab-inline"><span class="ab-mini">B</span>${vb}</span></dd>`;
    else dd = `<dd class="ab-rows"><span class="ab-row"><span class="ab-mini">A</span><span>${va}</span></span><span class="ab-row"><span class="ab-mini">B</span><span>${vb}</span></span></dd>`;
    return `<dt>${label}</dt>${dd}`;
  }).join('');
  $('sum-dataset').textContent = (lattice ? `${a.toFixed(3)} ${b.toFixed(3)} ${c.toFixed(3)} Å · ${gamma.toFixed(1)}°` : `${bins.join('×')}`)
    + (compare?.ready ? ' · vs B' : '');
  $('pipe-measured').textContent = `${(stats.valid / 1e6).toFixed(1)} M · ${pct(stats.fraction)} of the grid`;

  // Full details in the info popover.
  const items = [
    ['File', `${sourceName} (${mb(sourceSize)})`],
    ['Signal', meta.signal],
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
  if (compare?.ready) {
    const B = compare.meta, l = B.lattice;
    items.push(['Compare (B)', `${compare.name} (${mb(compare.size)})\n${[0, 1, 2].map((d) => B.shape[2 - d]).join(' × ')} bins, ${pct(B.stats.fraction)} measured`
      + (l ? `\na ${l.a.toFixed(4)}, b ${l.b.toFixed(4)}, c ${l.c.toFixed(4)} Å, γ ${l.gamma.toFixed(3)}°` : '')
      + (compare.mask ? `\nmask: ${pct((compare.mask.edge + compare.mask.outlier) / compare.mask.measured)} removed` : '')]);
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

/** The Dataset table's values (HTML) for one dataset. */
function datasetFacts(m, userMask) {
  const l = m.lattice, r = recipOf(l), bins = [0, 1, 2].map((d) => m.shape[2 - d]);
  const widths = m.dims.map((d) => fmt((d.edges[d.edges.length - 1] - d.edges[0]) / (d.edges.length - 1), 4));
  const same = (xs) => xs.every((x) => x === xs[0]);
  return {
    cell: l ? `${l.a.toFixed(3)} ${l.b.toFixed(3)} ${l.c.toFixed(3)} Å<br>${l.alpha.toFixed(2)}° ${l.beta.toFixed(2)}° ${l.gamma.toFixed(2)}° <span class="note">(${escapeHTML(l.source)})</span>` : 'not in file',
    recip: r ? `${r.len.map((x) => x.toFixed(4)).join(' ')} Å⁻¹<br>${r.angles.map((x) => x.toFixed(2)).join('° ')}°` : '—',
    grid: `${same(bins) ? `${bins[0]}³` : bins.join(' × ')}, Δ ${same(widths) ? widths[0] : widths.join(' ')}`,
    measured: pct(m.stats.fraction),
    masked: userMask ? pct((userMask.edge + userMask.outlier) / userMask.measured) : 'none',
  };
}

/** Processing pipeline, section badges and the legend above the views. */
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
  $('pipe-views').textContent = compare?.ready ? '3 slices, A | B split + 3-D (A)' : '3 slices + 3-D';
  $('sum-processing').textContent = `${mask ? `mask ${removed}` : 'no mask'} · ${sym ? symmetry.name : 'no symmetry'}`;
  $('sum-display').textContent = `${$('cmap').value} · ${$('vmin').value}–${$('vmax').value} · ${$('scale').dataset.value}`;
  $('legend-min').textContent = fmtValue(Number($('vmin').value) || 0);
  $('legend-max').textContent = fmtValue(Number($('vmax').value) || 0);
  $('legend-scale').textContent = $('scale').dataset.value;
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
  focusBtn.className = maxBtn.className = 'icon-btn';
  focusBtn.type = maxBtn.type = 'button';
  focusBtn.onclick = () => setLayout(layout === 'focus' && primary === key ? 'quad' : 'focus', key);
  maxBtn.onclick = () => setLayout(layout === 'single' && primary === key ? lastMulti : 'single', key);
  section.querySelector('.view-head').ondblclick = (e) => { if (!e.target.closest('button')) maxBtn.click(); };
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

function setupViewer() {
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
      fixed, x, y, step, lo, hi, version: 0, key, section: shell.section, zoom: null, drag: null, b: newLayer(),
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
    p.canvas.onpointerdown = (ev) => startBox(p, ev);
    p.canvas.onpointermove = (ev) => { hover(p, ev); moveBox(p, ev); };
    p.canvas.onpointerup = (ev) => endBox(p, ev);
    p.canvas.onpointercancel = () => { p.drag = null; p.canvas.classList.remove('panning'); redraw(); };
    p.canvas.onpointerleave = () => { p.hover.hidden = true; };
    p.canvas.ondblclick = () => { clearTimeout(p.clickTimer); setZoom(p, null); };
    p.zoomReset.onclick = () => setZoom(p, null);
    q('.save').onclick = () => savePNG(p);
    resized.observe(p.canvas);
    panels.push(p);
  }

  const shell = viewShell('3d', '<span class="badge">3D</span>', 'Isosurface');
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
    <input class="iso-level num" type="number" step="any" style="width: 76px" aria-label="Isosurface level">
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

  setLayout(layout, primary);
  paintAll();
}

// ---- Layout ------------------------------------------------------------------------

// Below this size the views stack in one column and the layouts do not apply.
const compactLayout = matchMedia('(max-width: 1000px), (max-height: 640px)');
compactLayout.addEventListener('change', () => { if (views.hk) setLayout(layout, primary); });

function setLayout(mode, key = primary) {
  if (!views[key]) key = 'hk';
  if (mode !== 'single') lastMulti = mode;
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
  }
  try { localStorage.setItem('nxv-layout', JSON.stringify({ mode: lastMulti, key })); } catch { /* storage unavailable */ }
}

// Redraw a plot whenever its canvas changes size (layout switches, window resizes).
let resizeQueued = false;
const resized = new ResizeObserver(() => {
  if (resizeQueued) return;
  resizeQueued = true;
  requestAnimationFrame(() => { resizeQueued = false; redraw(); });
});

// ---- Box zoom and clicks -------------------------------------------------------------

function localPoint(p, e) {
  const box = p.canvas.getBoundingClientRect();
  return [e.clientX - box.left, e.clientY - box.top];
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
    angles: $('angles').checked ? 'nominal' : 'measured', guides: $('guides').checked,
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
  const positive = [];
  for (const p of panels) {
    for (const layer of compare?.ready ? [p, p.b] : [p]) for (const v of layer.data?.values ?? []) if (v > 0) positive.push(v);
  }
  positive.sort((a, b) => a - b);
  const quantile = (q) => positive[Math.min(positive.length - 1, Math.floor(q * positive.length))];
  // Bragg peaks dominate the top percentiles, so a 97th-percentile ceiling
  // with the median as asinh softening keeps diffuse intensity visible.
  const vmax = positive.length ? sig(quantile(0.97)) : 1;
  $('vmin').value = $('scale').dataset.value === 'log' ? sig(vmax / 1000) : 0;
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
  for (const id of ['cmap-bar', 'legend-bar']) {
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

/** The colored image of one dataset's slice (a panel, or its B layer `p.b`), cached. */
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

  // Intensity images, clipped to the view parallelogram. When comparing, the
  // diagonal from its top-left to its bottom-right corner splits it: A below, B above.
  const corners = box.map(([u, v]) => project(u, v));
  const path = (pts) => {
    c.beginPath();
    pts.forEach(([a, b], i) => (i ? c.lineTo(a, b) : c.moveTo(a, b)));
    c.closePath();
  };
  path(corners);
  c.fillStyle = MISSING;
  c.fill();
  const shown = compareShown();
  if (shown === 'split' || shown === 'b') {
    // B's half is hatched parallel to the cut, so it stands out where it has no data.
    c.save();
    path(shown === 'split' ? corners.slice(1) : corners);
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
  const layers = shown === 'split' ? [[p, meta.dims, [corners[0], corners[1], corners[3]]], [p.b, compare.meta.dims, corners.slice(1)]]
    : shown === 'b' ? [[p.b, compare.meta.dims, corners]] : [[p, meta.dims, corners]];
  for (const [layer, dims, clip] of layers) {
    if (!layer.data) continue;
    const { rows, cols } = layer.data, ex = dims[p.x].edges, ey = dims[p.y].edges;
    const dx = (ex[ex.length - 1] - ex[0]) / cols, dy = (ey[ey.length - 1] - ey[0]) / rows;
    c.save();
    path(clip);
    c.clip();
    const [x0, y0] = project(ex[0], ey[0]);
    c.transform(sx * g.lx * dx, 0, sx * g.ly * g.cos * dy, -sy * g.ly * sin * dy, x0, y0);
    c.imageSmoothingEnabled = false;
    c.drawImage(layerImage(layer, s), 0, 0);
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

  // The cut along the diagonal (a white gap with dark edges) and dataset tags when comparing.
  if (shown === 'split') {
    c.save();
    c.lineCap = 'butt';
    c.shadowColor = 'rgba(18, 24, 33, 0.25)';
    c.shadowBlur = 3;
    for (const [color, width] of [['rgba(18, 24, 33, 0.75)', 4.5], ['#ffffff', 2.5]]) {
      c.strokeStyle = color;
      c.lineWidth = width;
      c.beginPath(); c.moveTo(...corners[3]); c.lineTo(...corners[1]); c.stroke();
      c.shadowColor = 'transparent';
    }
    c.restore();
  }
  if (shown) {
    const room = 0.45 * Math.hypot(corners[1][0] - corners[0][0], corners[1][1] - corners[0][1]);
    if (shown !== 'b') datasetTag(c, 'A', sourceName, corners[0], corners[1], corners[3], room, false);
    if (shown !== 'a') datasetTag(c, 'B', compare.name, corners[2], corners[3], corners[1], room, true);
  }

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
  const d = shown === 'b' && p.b.data ? p.b.data : p.data;
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
 * A dataset tag (letter and file name) inside the view corner `at`, offset
 * toward its neighbouring corners `along` (same row) and `side` (same column).
 * The tag grows up and right from a bottom corner, or down and left when `flip`.
 */
function datasetTag(c, letter, name, at, along, side, room, flip) {
  const unit = ([x, y]) => {
    const n = Math.hypot(x - at[0], y - at[1]) || 1;
    return [(x - at[0]) / n, (y - at[1]) / n];
  };
  const [ax, ay] = unit(along), [bx, by] = unit(side);
  const x = at[0] + 12 * ax + 8 * bx, y = at[1] + 12 * ay + 8 * by;
  c.save();
  c.font = `600 12px ${SANS}`;
  const h = 26, badge = 20, maxText = room - badge - 16;
  // A long name keeps the part that differs from the other dataset's name.
  const fits = (t) => c.measureText(t).width <= maxText;
  const text = maxText > 24 ? shortenName(name, fits, name === sourceName ? compare?.name : sourceName) : '';
  const tw = text && fits(text) ? c.measureText(text).width : 0;
  const w = 3 + badge + (tw ? 7 + tw + 9 : 3), left = flip ? x - w : x, top = flip ? y : y - h;
  c.fillStyle = 'rgba(255, 255, 255, 0.9)';
  c.strokeStyle = 'rgba(18, 24, 33, 0.18)';
  c.lineWidth = 1;
  c.beginPath(); c.roundRect(left, top, w, h, 6); c.fill(); c.stroke();
  c.fillStyle = INK;
  c.beginPath(); c.roundRect(left + 3, top + 3, badge, badge, 5); c.fill();
  c.fillStyle = '#ffffff';
  c.textAlign = 'center';
  c.textBaseline = 'middle';
  c.font = `700 12.5px ${SANS}`;
  c.fillText(letter, left + 3 + badge / 2, top + h / 2 + 0.5);
  if (tw) {
    c.font = `600 12px ${SANS}`;
    c.fillStyle = INK2;
    c.textAlign = 'left';
    c.fillText(text, left + 3 + badge + 7, top + h / 2 + 0.5);
  }
  c.restore();
}

/** Redraw one panel, with the zoom box being dragged (as its u-v parallelogram). */
function drawPanel(p) {
  if (!settings) return;
  draw(p, p.canvas, p.canvas.clientWidth, p.canvas.clientHeight, devicePixelRatio || 1);
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
    p.angle.title = `Angle between the ${meta.dims[p.x].label} and ${meta.dims[p.y].label} axes (reciprocal lattice, ${settings.angles} cell angles)`;
    drawPanel(p);
  }
  paintColorbar();
  updateStates();
  update3D();
}

// ---- Interaction -----------------------------------------------------------------

function pointAt(p, e) {
  if (!p.data || !p.inverse) return null;
  const box = p.canvas.getBoundingClientRect();
  const [u, v] = p.inverse(e.clientX - box.left, e.clientY - box.top);
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
    // Both datasets, the one under the cursor first.
    const a = `A ${read(p, meta.dims) ?? 'no data'}`, b = `B ${read(p.b, compare.meta.dims) ?? 'no data'}`;
    const [u0, u1, v0, v1] = p.view;
    const inB = shown === 'b' || (shown === 'split' && (pt[0] - u0) / (u1 - u0) + (pt[1] - v0) / (v1 - v0) > 1);
    text = inB ? `${b} · ${a}` : `${a} · ${b}`;
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
  const out = document.createElement('canvas');
  draw(p, out, Math.max(p.canvas.clientWidth, 480), Math.max(p.canvas.clientHeight, 360) + 24, 3, true);
  const X = meta.dims[p.x], Y = meta.dims[p.y], F = meta.dims[p.fixed];
  const names = { a: stem(), b: stemOf(compare?.name ?? ''), split: `${stem()}_vs_${stemOf(compare?.name ?? '')}` };
  out.toBlob((blob) => download(blob, `${names[compareShown() ?? 'a']}_${X.label}${Y.label}_${F.label}=${p.data.center}.png`));
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
  if (compare?.ready) {
    Object.assign(compare, compareMaps());
    showCompare();
  }
  showSymmetry();
  describe();
  panels.forEach(request);
  requestIso();
}

function showSymmetry() {
  const { ops } = symmetry, statusEl = $('sym-status');
  statusEl.className = 'note';
  if (ops.length === 1) {
    statusEl.textContent = 'No symmetry averaging: each voxel is used as measured.';
  } else {
    let text = `${ops.length} operations; equivalent voxels are pooled with equal weight.`;
    const hkl = meta.dims.every((d) => d.basis);
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
      text += ' Axes have no HKL basis, so operations act on the display axes directly.';
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
  status('busy', 'building mask');
  $('mask-state').textContent = 'working…';
  $('mask-state').className = 'state busy';
  $('mask-status').className = 'note';
  $('mask-status').textContent = radius || k ? 'Building mask…' : 'Clearing mask…';
  worker.postMessage({ type: 'mask', id: ++requestId, radius, k, maps: symmetry.maps, symmetry: symmetry.name });
  if (compare?.ready) sendCompareMask(radius, k);
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
  const statusEl = $('export-status'), state = $('export-state');
  $('export-progress').hidden = true;
  let plan;
  try {
    plan = exportPlan(meta.dims, meta.lattice);
  } catch (err) {
    $('export-run').disabled = true;
    state.textContent = 'unavailable';
    state.className = 'state';
    statusEl.className = 'note error';
    statusEl.textContent = err.message;
    return;
  }
  $('export-run').disabled = false;
  state.textContent = '3D-ΔPDF';
  state.className = 'state';
  const sym = symmetry.ops.length > 1, voxels = plan.shape.reduce((a, b) => a * b);
  const parts = [
    `${plan.shape.join(' × ')} (H × K × L)${plan.padded ? ', padded to be symmetric about 0' : ''}`,
    sym ? `${symmetry.name} averaged` : 'not symmetrized',
    mask ? `mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : 'no mask',
  ];
  if (compare?.ready) parts.push('dataset A');
  const warnings = [];
  if (!sym) warnings.push('NEBULA3D expects a symmetrized volume: choose a Laue class first.');
  if (voxels > 80e6) warnings.push(`NEBULA3D in the browser handles up to about 80 M voxels; this is ${Math.round(voxels / 1e6)} M, so use its desktop app.`);
  statusEl.className = warnings.length ? 'note warn' : 'note';
  statusEl.textContent = `${parts.join(' · ')}.${warnings.length ? ` ${warnings.join(' ')}` : ''}`;
}

function showExportProgress(label, fraction) {
  $('export-progress').hidden = false;
  $('export-bar').style.width = `${Math.round(100 * clamp(fraction, 0, 1))}%`;
  $('export-status').className = 'note';
  $('export-status').textContent = `${label}… ${Math.round(100 * fraction)}%`;
}

function runExport() {
  if (!panels.length || exportJob) return;
  let plan, maps;
  try {
    plan = exportPlan(meta.dims, meta.lattice);
    maps = symmetry.ops.length > 1 ? indexMaps(symmetry.ops, plan.paddedDims) : [IDENTITY_MAP];
  } catch (err) {
    $('export-status').className = 'note error';
    $('export-status').textContent = err.message;
    return;
  }
  const sym = symmetry.ops.length > 1;
  exportJob = { name: `${stem()}_${sym ? `sym${symmetry.name.replace(/\//g, '')}` : 'unsym'}.nxs` };
  $('export-run').disabled = true;
  $('export-state').textContent = 'working…';
  $('export-state').className = 'state busy';
  status('busy', 'exporting');
  showExportProgress('Symmetrizing', 0);
  const { order, lo, size, shape, centers, ub, padded } = plan;
  worker.postMessage({
    type: 'export', id: ++requestId, plan: { order, lo, size, shape, centers, ub }, maps,
    attrs: {
      source_file: sourceName,
      symmetry: sym ? symmetry.name : 'none',
      symmetry_ops: symmetry.ops.map(formatOp).join('; '),
      mask: mask ? `coverage-edge erosion ${mask.radius}, outlier cut ${mask.k} sigma, ${pct((mask.edge + mask.outlier) / mask.measured)} of measured voxels removed` : 'none',
      ub_source: meta.lattice.source,
      padded: padded ? 'yes' : 'no',
      created_by: 'NeXus Viewer, https://drthyang.github.io/neutron-nexus-viewer/',
      created: new Date().toISOString(),
    },
  });
}

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

function update3D() {
  if (!view3d || !settings) return;
  const cell = meta.lattice && (settings.angles === 'nominal' ? nominalCell(meta.lattice) : meta.lattice);
  const { T } = cartesianBasis(meta.dims, cell);
  view3d.setFrame(T, meta.dims.map((d) => viewRange(d, settings.limit)), meta.dims.map((d) => d.label));
  const slices = panels.filter((p) => p.data).map((p) => {
    const X = meta.dims[p.x], Y = meta.dims[p.y], u = viewRange(X, settings.limit), v = viewRange(Y, settings.limit);
    return {
      fixed: p.fixed, x: p.x, y: p.y, center: p.data.center, ...sliceTexture(p, settings, u, v), u, v,
      ex: [X.edges[0], X.edges[X.edges.length - 1]], ey: [Y.edges[0], Y.edges[Y.edges.length - 1]],
    };
  });
  view3d.setSlices(slices, $('iso-slices').checked);
}

/**
 * The texture of a slice plane in the 3-D view, on A's pixel grid. When
 * comparing, it matches the 2-D view: split along the diagonal of the visible
 * range (u, v), or B alone.
 */
function sliceTexture(p, s, u, v) {
  const shown = compareShown(), a = layerImage(p, s);
  if (!shown || shown === 'a') return { image: a, key: p.imageKey };
  const b = p.b.data ? layerImage(p.b, s) : null;
  const key = `${p.imageKey}|${b ? p.b.imageKey : '-'}|${shown}|${u}|${v}`;
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
  if (shown === 'split') {
    c.save();
    clip([px(u[0], v[0]), px(u[1], v[0]), px(u[0], v[1])]);
    c.drawImage(a, 0, 0);
    c.restore();
  }
  if (b) {
    const bx = compare.meta.dims[p.x].edges, by = compare.meta.dims[p.y].edges, { rows: rb, cols: cb } = p.b.data;
    c.save();
    if (shown === 'split') clip([px(u[1], v[0]), px(u[1], v[1]), px(u[0], v[1])]);
    c.transform((bx[bx.length - 1] - bx[0]) / cb / dx, 0, 0, (by[by.length - 1] - by[0]) / rb / dy, (bx[0] - ex[0]) / dx, (by[0] - ey[0]) / dy);
    c.imageSmoothingEnabled = false;
    c.drawImage(b, 0, 0);
    c.restore();
  }
  if (shown === 'split') {
    c.strokeStyle = '#ffffff';
    c.lineWidth = Math.max(1, Math.max(cols, rows) / 120);
    c.beginPath(); c.moveTo(...px(u[0], v[1])); c.lineTo(...px(u[1], v[0])); c.stroke();
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
  pendingCompare = DEMO.compare;
  pendingProcessing = { sym: DEMO.sym, mask: DEMO.mask };
  openURL(DEMO.url);
}

// ---- Comparing two datasets ---------------------------------------------------------------

/** Open a second file (B) in its own worker; its slices follow A's positions and processing. */
function openCompare(file) {
  if (!panels.length) return; // viewer not ready yet
  closeCompare();
  const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  const current = compare = {
    worker: w, name: file.name, size: file.size, meta: null, ready: false, autoscaled: false,
    maps: [IDENTITY_MAP], mapsNote: '', mask: null, maskNote: '', maskBusy: false,
  };
  w.onmessage = ({ data }) => { if (compare === current) compareHandlers[data.type]?.(data); };
  w.onerror = (e) => {
    e.preventDefault();
    if (compare !== current) return;
    closeCompare();
    error(`Dataset B: ${e.message || 'the HDF5 reader could not start.'}`);
  };
  w.postMessage({ type: 'open', file });
  showCompare('Opening', 0);
}

async function openCompareURL(url) {
  if (!panels.length) return;
  closeCompare();
  const current = compare = { worker: null, name: urlName(url), ready: false };
  showCompare('Downloading', 0);
  try {
    const file = await fetchFile(url, (got, total) => { if (compare === current) showCompare(`Downloading ${downloaded(got, total)}`, total ? got / total : 0); });
    if (compare === current) openCompare(file);
  } catch (err) {
    if (compare !== current) return;
    closeCompare();
    error(`Could not download ${url}: ${err.message}. The server must allow cross-origin requests.`);
  }
}

function closeCompare() {
  compare?.worker?.terminate();
  compare = null;
  for (const p of panels) {
    p.b = newLayer();
    showCaption(p);
  }
  showCompare();
  if (meta && panels.length) {
    describe();
    redraw();
  }
}

/** The symmetry's index maps on B's grid, or the identity (with a note) when they do not fit it. */
function compareMaps() {
  if (symmetry.ops.length === 1) return { maps: [IDENTITY_MAP], mapsNote: '' };
  try {
    return { maps: indexMaps(symmetry.ops, compare.meta.dims), mapsNote: '' };
  } catch (err) {
    return { maps: [IDENTITY_MAP], mapsNote: `Symmetry not applied to B: ${err.message}` };
  }
}

function sendCompareMask(radius, k) {
  const kb = compare.maps.length >= 3 ? k : 0;
  compare.maskNote = k && !kb ? 'Outlier cut not applied to B: it needs a symmetry of at least 3 operations that fits B\'s grid.' : '';
  compare.maskBusy = true;
  compare.maskRequested = `${radius},${k}`;
  compare.worker.postMessage({ type: 'mask', id: ++requestId, radius, k: kb, maps: compare.maps, symmetry: symmetry.name });
  showCompare('Building mask', 0);
}

const compareHandlers = {
  progress: ({ label, fraction }) => showCompare(label, fraction),
  meta: ({ info }) => { compare.meta = info; },
  ready: ({ stats, seconds }) => {
    Object.assign(compare.meta, { stats, seconds });
    compare.ready = true;
    Object.assign(compare, compareMaps());
    for (const p of panels) p.b = newLayer();
    // B gets A's mask before its first slices.
    if (mask) sendCompareMask(mask.radius, mask.k);
    else panels.forEach((p) => request(p, 'b'));
    if (compare.maskBusy) showCompare('Building mask', 0);
    else showCompare();
    describe();
    redraw();
  },
  slice: (msg) => {
    const p = panels.find((q) => q.fixed === msg.fixed);
    if (!p) return;
    Object.assign(p.b, { busy: false, data: msg, error: '', version: p.b.version + 1 });
    showCaption(p);
    if (p.b.wanted) sendB(p);
    // Once B's first slices are in, the automatic range covers both datasets.
    if (!compare.autoscaled && panels.every((q) => q.b.data || q.b.error)) {
      compare.autoscaled = true;
      if (rangeIsAuto) autoRange();
    }
    redraw();
  },
  error: ({ fixed, message }) => {
    const p = panels.find((q) => q.fixed === fixed);
    if (!p) {
      closeCompare();
      error(`Dataset B: ${message}`);
      return;
    }
    Object.assign(p.b, { busy: false, data: null, error: message, version: p.b.version + 1 });
    showCaption(p);
    if (p.b.wanted) sendB(p);
    redraw();
  },
  'progress-mask': ({ label, fraction }) => showCompare(label, fraction),
  mask: ({ stats, radius, k }) => {
    compare.maskBusy = false;
    compare.mask = stats ? { ...stats, radius, k } : null;
    showCompare();
    describe();
    panels.forEach((p) => request(p, 'b'));
  },
  'mask-error': ({ message }) => {
    compare.maskBusy = false;
    compare.maskNote = `Mask for B failed: ${message}`;
    showCompare();
    panels.forEach((p) => request(p, 'b'));
  },
};

/**
 * The Compare card, the B button in the top bar and the A / Split / B control.
 * `step` and `fraction` describe work in progress on B (loading or masking).
 */
function showCompare(step = '', fraction = null) {
  const ready = !!compare?.ready;
  document.body.classList.toggle('comparing', ready);
  // Loaded, A and B share the Dataset table; the Compare card only opens and loads B.
  $('compare').hidden = ready;
  $('data-files').hidden = !ready;
  $('compare-open').hidden = !!compare;
  $('compare-file').hidden = !compare || ready;
  $('compare-progress').hidden = !compare || ready || fraction === null;
  if (sourceName) setFileName($('dataset-name'), sourceName, ready ? compare.name : null);
  if (compare && !ready) {
    setFileName($('compare-name'), compare.name, sourceName);
    $('compare-size').textContent = compare.size ? mb(compare.size) : '';
    $('compare-bar').style.width = `${Math.round(100 * clamp(fraction ?? 0, 0, 1))}%`;
  }
  const stepText = step && `${step}${fraction > 0 ? ` · ${Math.round(100 * fraction)}%` : '…'}`;
  $('compare-view').hidden = $('compare-button').hidden = $('file-sep').hidden = $('dataset-tag').hidden = !ready;
  $('compare-state').textContent = compare ? 'loading' : 'off';
  $('compare-state').className = `state${compare ? ' busy' : ''}`;
  const title3d = views['3d']?.section.querySelector('.view-title');
  if (title3d) title3d.textContent = ready ? 'Isosurface · A' : 'Isosurface';
  const note = $('compare-note');
  note.hidden = true;
  if (!compare) {
    $('compare-status').textContent = 'Split every slice along its diagonal: this dataset (A) below, a second one (B) above.';
    return;
  }
  if (!ready) {
    $('compare-status').textContent = stepText || 'Opening…';
    return;
  }
  setFileName($('file-a-name'), sourceName, compare.name);
  $('file-a-size').textContent = mb(sourceSize);
  setFileName($('file-b-name'), compare.name, sourceName);
  $('file-b-size').textContent = compare.size ? mb(compare.size) : '';
  setFileName($('compare-button-name'), compare.name, sourceName);
  $('compare-button').title = `Dataset B: ${compare.name} (click to replace)`;
  // Work in progress on B, or warnings, under the Dataset table.
  const labels = (dims) => dims.map((d) => d.label).join(' ');
  const warnings = [];
  if (labels(compare.meta.dims) !== labels(meta.dims)) warnings.push(`B's axes (${labels(compare.meta.dims)}) differ from A's (${labels(meta.dims)}); B is drawn on A's axes.`);
  if (compare.mapsNote) warnings.push(compare.mapsNote);
  if (compare.maskNote) warnings.push(compare.maskNote);
  if (stepText || warnings.length) {
    note.hidden = false;
    note.className = stepText ? 'note' : 'note warn';
    note.textContent = stepText ? `B: ${stepText}` : warnings.join(' ');
  }
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
  if (typeof saved.angles === 'boolean') $('angles').checked = saved.angles;
  if (typeof saved.guides === 'boolean') $('guides').checked = saved.guides;
  setClickMode(['zoom', 'move'].includes(saved.clickMode) ? saved.clickMode : 'navigate');
  try {
    const l = JSON.parse(localStorage.getItem('nxv-layout'));
    if (['quad', 'focus'].includes(l?.mode)) lastMulti = layout = l.mode;
    if (['hk', 'hl', 'kl', '3d'].includes(l?.key)) primary = l.key;
  } catch { /* storage unavailable: keep the defaults */ }
  paintColorbar();
}

function persist() {
  const values = { cmap: $('cmap').value, scale: $('scale').dataset.value, angles: $('angles').checked, guides: $('guides').checked, clickMode };
  try { localStorage.setItem('nxv-settings', JSON.stringify(values)); } catch { /* storage unavailable */ }
}

restore();
show('intro');
status('', 'no file');
$('open').onclick = $('open-intro').onclick = () => $('file').click();
$('open-example').onclick = openDemo;
$('file').onchange = () => {
  if ($('file').files[0]) {
    pendingCompare = pendingProcessing = null;
    openFile($('file').files[0]);
  }
  $('file').value = '';
};
$('compare-open').onclick = $('compare-replace').onclick = $('compare-button').onclick = () => $('file-b').click();
$('file-b').onchange = () => { if ($('file-b').files[0]) openCompare($('file-b').files[0]); $('file-b').value = ''; };
$('compare-close').onclick = $('compare-cancel').onclick = closeCompare;
segmented($('compare-view'), (value) => {
  compareView = value;
  redraw();
});
$('error-close').onclick = () => error('');
for (const id of ['cmap', 'vmin', 'vmax', 'soft', 'limit', 'angles', 'guides']) $(id).addEventListener('input', redraw);
for (const id of ['cmap', 'angles', 'guides']) $(id).addEventListener('change', persist);
for (const id of ['vmin', 'vmax', 'soft']) $(id).addEventListener('input', () => { rangeIsAuto = false; });
$('angles').addEventListener('change', () => { if (meta) describe(); });
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
for (const li of document.querySelectorAll('.pipeline li[data-target]')) {
  li.onclick = () => {
    setPanel(true);
    const sec = li.closest('.psec');
    if (sec.classList.contains('collapsed')) sec.querySelector('.psec-head').click();
    const target = $(li.dataset.target);
    target.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    target.classList.remove('flash');
    void target.offsetWidth;
    target.classList.add('flash');
  };
}
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
$('export-run').onclick = runExport;

// Dropping a file opens it, or opens it as dataset B over the Compare card or the B button.
const dropsOnB = (e) => panels.length > 0 && !!e.target.closest?.('#compare, #compare-button, #data-files');
document.addEventListener('dragover', (e) => {
  e.preventDefault();
  document.body.classList.add('dragging');
  document.body.classList.toggle('drop-b', dropsOnB(e));
});
document.addEventListener('dragleave', (e) => { if (!e.relatedTarget) document.body.classList.remove('dragging', 'drop-b'); });
document.addEventListener('drop', (e) => {
  e.preventDefault();
  document.body.classList.remove('dragging', 'drop-b');
  const file = e.dataTransfer?.files?.[0];
  if (!file) return;
  if (dropsOnB(e)) openCompare(file);
  else {
    pendingCompare = pendingProcessing = null;
    openFile(file);
  }
});

// ?url= opens a remote file, ?compare= a second one as dataset B, and ?sym=
// (a Laue class, such as 6/mmm) and ?mask= (erosion radius, optionally ",k"
// for the outlier cut) process them. ?demo opens the example.
const params = new URLSearchParams(location.search), remote = params.get('url');
if (params.has('demo')) openDemo();
else if (remote) {
  pendingCompare = params.get('compare');
  pendingProcessing = { sym: params.get('sym'), mask: params.get('mask') };
  openURL(remote);
}
