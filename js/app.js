import { COLORMAPS } from './colormaps.js';
import { cartesianBasis, nominalCell, planeGeometry } from './nexus.js';
import { IDENTITY_MAP } from './slab.js';
import { closeGroup, formatOp, indexMaps, metricChange, parseOps, PRESETS } from './symmetry.js';

const $ = (id) => document.getElementById(id);
// [fixed, x, y] display dimensions: HK, HL and KL planes for Mantid HKL data.
const LAYOUT = [[2, 0, 1], [1, 0, 2], [0, 1, 2]];
const LUTS = Object.fromEntries(Object.entries(COLORMAPS).map(([name, hex]) =>
  [name, Uint8Array.from(hex.match(/../g), (h) => parseInt(h, 16))]));
const SAVED = ['cmap', 'scale', 'angles', 'guides'];
const ICON_ENLARGE = '<svg width="13" height="13" viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M8.5 1.5h4v4M5.5 12.5h-4v-4M12.5 1.5 8 6M1.5 12.5 6 8"/></svg>';
const ICON_RESTORE = '<svg width="13" height="13" viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M12.5 5.5h-4v-4M1.5 8.5h4v4M8.5 5.5 13 1M5.5 8.5 1 13"/></svg>';
const NO_SYMMETRY = { name: '1', ops: [[1, 0, 0, 0, 1, 0, 0, 0, 1]], maps: [IDENTITY_MAP] };
// Canvas colors, matching the page tokens.
const INK = '#141a22', INK2 = '#4a5566', AXIS = '#9aa5b3', MISSING = '#e7ebf0';
const MONO = 'ui-monospace, "SF Mono", Menlo, Consolas, monospace';
const SANS = '-apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, sans-serif';

let worker = null, meta = null, panels = [], settings = null, sourceName = '', sourceSize = 0, autoscaled = false;
let requestId = 0, symmetry = NO_SYMMETRY, mask = null;
// Slice layout: the fixed dimension of the enlarged plane, or null for three across.
let focus = 2;
// 3-D view state: the lazily loaded View3D, its panel elements and the isosurface request queue.
let view3d = null, iso = null;

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

function niceTicks(lo, hi, count = 5) {
  const raw = (hi - lo) / count;
  if (!(raw > 0)) return [];
  const p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw * 0.999);
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

function engine(state, text) {
  $('engine-dot').className = `engine-dot ${state}`;
  $('engine-text').textContent = text;
}

// ---- Files and the worker ---------------------------------------------------

function show(stage) {
  for (const id of ['intro', 'loading', 'viewer']) $(id).hidden = id !== stage;
}

function progress(label, fraction) {
  $('progress-label').textContent = label;
  $('progress').style.width = `${Math.round(100 * clamp(fraction, 0, 1))}%`;
}

function error(message) {
  $('error').textContent = message;
  $('error').hidden = !message;
}

function fail(message) {
  error(message);
  engine(meta ? 'ok' : '', meta ? 'in-browser engine' : 'no file open');
  if (!meta || $('viewer').hidden) show('intro');
}

function openFile(file) {
  worker?.terminate();
  meta = null;
  panels = [];
  autoscaled = false;
  symmetry = NO_SYMMETRY;
  mask = null;
  iso = null;
  sourceName = file.name;
  sourceSize = file.size;
  $('slices').replaceChildren();
  $('view3d').replaceChildren();
  error('');
  show('loading');
  engine('busy', 'reading file');
  $('dataset-name').textContent = file.name;
  progress(`Opening ${file.name} (${mb(file.size)})`, 0);
  document.title = `${file.name} · NeXus Slice Viewer`;
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = ({ data }) => handlers[data.type]?.(data);
  worker.onerror = (e) => {
    e.preventDefault();
    fail(e.message || 'The HDF5 reader could not start. It is loaded from cdn.jsdelivr.net; check the network connection.');
  };
  worker.postMessage({ type: 'open', file });
}

async function openURL(url) {
  show('loading');
  engine('busy', 'downloading');
  progress(`Downloading ${url}`, 0);
  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const total = Number(res.headers.get('content-length')) || 0;
    const reader = res.body.getReader(), parts = [];
    for (let got = 0; ;) {
      const { done, value } = await reader.read();
      if (done) break;
      parts.push(value);
      got += value.length;
      progress(`Downloading ${mb(got)}${total ? ` of ${mb(total)}` : ''}`, total ? got / total : 0);
    }
    const name = decodeURIComponent(new URL(url, location.href).pathname.split('/').pop()) || 'remote.nxs';
    openFile(new File(parts, name));
  } catch (err) {
    fail(`Could not download ${url}: ${err.message}. The server must allow cross-origin requests.`);
  }
}

const handlers = {
  progress: ({ label, fraction }) => progress(`${label}… ${Math.round(100 * fraction)}%`, fraction),
  meta: ({ info }) => { meta = info; },
  ready: ({ stats, seconds }) => {
    meta.stats = stats;
    meta.seconds = seconds;
    engine('ok', 'in-browser engine');
    setupViewer();
    show('viewer');
    redraw();
    panels.forEach(request);
    setup3D();
  },
  slice: (msg) => {
    const p = panels.find((q) => q.fixed === msg.fixed);
    if (!p) return;
    p.busy = false;
    p.data = msg;
    p.version++;
    p.caption.textContent = `${meta.dims[p.fixed].label} ∈ [${fmt(msg.slab[0])}, ${fmt(msg.slab[1])}] · ${msg.bins} bin${msg.bins === 1 ? '' : 's'}`
      + ` · ${pct(msg.coverage)} measured${symmetryNote(msg)}${msg.removed ? ' · removed voxels only' : ''}`;
    p.tag.textContent = `${meta.dims[p.fixed].label} = ${fmt(msg.center)}`;
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
    iso.caption.textContent = `${msg.shape.join(' × ')} grid (${msg.factor}× binned) · level ${fmtValue(msg.level)} · `
      + `${triangles.toLocaleString()} triangles${msg.symmetry !== '1' ? ` · ${msg.symmetry}` : ''}${mask ? ' · masked' : ''}${msg.note ? ` · ${msg.note}` : ''}`;
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
    engine('ok', 'in-browser engine');
    $('mask-apply').disabled = false;
    mask = stats ? { ...stats, radius, k, group } : null;
    $('mask-clear').disabled = $('mask-download').disabled = $('mask-removed').disabled = !mask;
    if (!mask) $('mask-removed').checked = false;
    showMask(seconds);
    panels.forEach(request);
    requestIso();
    describe();
  },
  'mask-error': ({ message }) => {
    engine('ok', 'in-browser engine');
    $('mask-apply').disabled = false;
    $('mask-status').className = 'note error';
    $('mask-status').textContent = message;
  },
  'mask-file': ({ blob }) => download(blob, `${stem()}_mask.npy.gz`),
};

function request(p) {
  const center = Number(p.center.value), thickness = Number(p.width.value);
  if (!Number.isFinite(center) || !Number.isFinite(thickness) || thickness <= 0) {
    p.caption.textContent = 'Enter a finite center and a positive thickness.';
    return;
  }
  p.wanted = { center, thickness };
  if (!p.busy) send(p);
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

const symmetryNote = (data) => (data.order > 1 ? ` · ${data.symmetry} (${data.order} ops)` : '');
const stem = () => sourceName.replace(/\.[^.]+$/, '');

// ---- Viewer setup -------------------------------------------------------------

function describe() {
  const { dims, shape, lattice, stats } = meta;
  const bins = [0, 1, 2].map((d) => shape[2 - d]);
  $('ident').innerHTML = dims.map((d, i) => `<span class="hue-${i}">${escapeHTML(d.label)}</span>`).join('<span class="dot">·</span>');
  $('ident-eyebrow').textContent = `${bins.join(' × ')} bins`;
  $('ident-desc').textContent = `${sourceName} · ${mb(sourceSize)} · ${pct(stats.fraction)} of voxels measured · read in ${meta.seconds.toFixed(1)} s`;
  const items = [
    ['Source', meta.signal],
    ['Axes', dims.map((d) => `${d.longName} ${fmt(d.edges[0], 2)}…${fmt(d.edges[d.edges.length - 1], 2)}`).join('  ')],
    ['Grid', `${bins.join(' × ')} · Δ ${dims.map((d) => fmt((d.edges[d.edges.length - 1] - d.edges[0]) / (d.edges.length - 1), 4)).join(', ')}`],
  ];
  if (lattice) {
    const { a, b, c, alpha, beta, gamma } = lattice;
    items.push([`Lattice (${lattice.source})`, `${a.toFixed(4)} ${b.toFixed(4)} ${c.toFixed(4)} Å · ${alpha.toFixed(2)} ${beta.toFixed(2)} ${gamma.toFixed(2)}°`]);
  }
  items.push(['Symmetry', symmetry.ops.length > 1 ? `${symmetry.name} · ${symmetry.ops.length} ops` : 'none']);
  items.push(['Mask', mask ? `${pct((mask.edge + mask.outlier) / mask.measured)} removed` : 'none']);
  $('meta').innerHTML = items.map(([k, v]) => `<div class="meta-item"><span class="meta-key">${escapeHTML(k)}</span><span class="meta-val" title="${escapeHTML(v)}">${escapeHTML(v)}</span></div>`).join('');
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
    const units = F.units ? ` ${escapeHTML(F.units)}` : '';
    const section = document.createElement('section');
    section.className = 'card';
    section.innerHTML = `
      <div class="card-head">
        <span class="card-titlegroup"><span class="badge hue-${fixed}">${escapeHTML(F.label)}</span>
          <span class="panel-title">${escapeHTML(X.label)} – ${escapeHTML(Y.label)}</span><span class="card-sub"></span></span>
        <span class="head-actions">
          <button type="button" class="btn btn-ghost btn-xs zoom-reset" hidden title="Back to the full view (or double-click the plot)">Reset zoom</button>
          <button type="button" class="btn btn-ghost btn-xs save-btn">Save PNG</button>
          <button type="button" class="icon-btn focus-btn"></button>
        </span>
      </div>
      <canvas class="plot" role="img"></canvas>
      <div class="cut cut--${fixed}">
        <div class="field">
          <div class="field-row"><span class="field-label">Center</span>
            <span class="readout-edit">${escapeHTML(F.label)} =<input class="center" type="number" step="any">${units}</span></div>
          <input class="slider" type="range">
        </div>
        <div class="field">
          <div class="field-row"><span class="field-label" title="Full slab thickness">Thickness</span>
            <span class="readout-edit"><input class="width" type="number" min="0" step="any">${units}</span></div>
          <input class="wslider" type="range">
        </div>
      </div>
      <div class="card-foot caption-foot"><span class="caption" aria-live="polite"></span><span class="hover">&nbsp;</span></div>`;
    $('slices').append(section);
    const q = (sel) => section.querySelector(sel);
    const p = {
      fixed, x, y, step, lo, hi, version: 0, section, zoom: null, drag: null, zoomReset: q('.zoom-reset'), focusBtn: q('.focus-btn'),
      slider: q('.slider'), center: q('.center'), wslider: q('.wslider'), width: q('.width'),
      caption: q('.caption'), tag: q('.card-sub'), canvas: q('canvas'), hover: q('.hover'),
    };
    p.canvas.setAttribute('aria-label', `${X.label}–${Y.label} intensity slice`);
    Object.assign(p.slider, { min: roundTo(lo, step), max: roundTo(hi, step), step, value: center });
    Object.assign(p.wslider, { min: step, max: roundTo(maxWidth, step), step, value: 3 * step });
    p.center.value = center;
    p.width.value = roundTo(3 * step, step);
    p.slider.oninput = () => { p.center.value = p.slider.value; request(p); };
    p.center.oninput = () => { p.slider.value = p.center.value; paint(p.slider); request(p); };
    p.wslider.oninput = () => { p.width.value = p.wslider.value; request(p); };
    p.width.oninput = () => { p.wslider.value = p.width.value; paint(p.wslider); request(p); };
    p.canvas.onpointerdown = (ev) => startBox(p, ev);
    p.canvas.onpointermove = (ev) => { hover(p, ev); moveBox(p, ev); };
    p.canvas.onpointerup = (ev) => endBox(p, ev);
    p.canvas.onpointercancel = () => { p.drag = null; redraw(); };
    p.canvas.onpointerleave = () => { p.hover.textContent = '\u00a0'; };
    p.canvas.ondblclick = () => { clearTimeout(p.clickTimer); setZoom(p, null); };
    p.zoomReset.onclick = () => setZoom(p, null);
    p.focusBtn.onclick = () => setFocus(focus === fixed ? null : fixed);
    q('.save-btn').onclick = () => savePNG(p);
    panels.push(p);
  }
  paintAll();
  setFocus(focus);
}

// ---- Layout and box zoom ---------------------------------------------------------

function setFocus(fixed) {
  focus = fixed;
  $('slices').classList.toggle('focus', fixed !== null);
  for (const p of panels) {
    const big = p.fixed === fixed;
    p.section.classList.toggle('is-focus', big);
    p.focusBtn.innerHTML = big ? ICON_RESTORE : ICON_ENLARGE;
    p.focusBtn.title = big ? 'Show the three slices side by side' : 'Enlarge this slice';
  }
  try { localStorage.setItem('nxv-focus', JSON.stringify(fixed)); } catch { /* storage unavailable */ }
  redraw();
}

function localPoint(p, e) {
  const box = p.canvas.getBoundingClientRect();
  return [e.clientX - box.left, e.clientY - box.top];
}

function startBox(p, e) {
  if (e.button !== 0 || !p.inverse) return;
  const [x, y] = localPoint(p, e);
  p.drag = { x0: x, y0: y, x1: x, y1: y, moved: false };
  try { p.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic or already released pointer */ }
}

function moveBox(p, e) {
  if (!p.drag) return;
  [p.drag.x1, p.drag.y1] = localPoint(p, e);
  p.drag.moved ||= Math.hypot(p.drag.x1 - p.drag.x0, p.drag.y1 - p.drag.y0) > 5;
  if (p.drag.moved) drawPanel(p);
}

function endBox(p, e) {
  const d = p.drag;
  p.drag = null;
  if (!d) return;
  if (!d.moved) {
    // Wait briefly so a double-click (reset zoom) does not also move the other slices.
    const { clientX, clientY } = e;
    clearTimeout(p.clickTimer);
    p.clickTimer = setTimeout(() => navigate(p, { clientX, clientY }), 250);
    return;
  }
  const box = boxWindow(p, d);
  if (box) setZoom(p, box);
  else drawPanel(p);
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
  return u[1] - u[0] >= wx && v[1] - v[0] >= wy ? { u, v } : null;
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
  for (const p of panels) for (const v of p.data?.values ?? []) if (v > 0) positive.push(v);
  positive.sort((a, b) => a - b);
  const quantile = (q) => positive[Math.min(positive.length - 1, Math.floor(q * positive.length))];
  // Bragg peaks dominate the top percentiles, so a 97th-percentile ceiling
  // with the median as asinh softening keeps diffuse intensity visible.
  const vmax = positive.length ? sig(quantile(0.97)) : 1;
  $('vmin').value = $('scale').dataset.value === 'log' ? sig(vmax / 1000) : 0;
  $('vmax').value = vmax;
  $('soft').value = positive.length ? sig(quantile(0.5), 1) : sig(vmax / 20);
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
  const c = $('cmap-bar').getContext('2d'), lut = LUTS[$('cmap').value];
  const img = c.createImageData(256, 1);
  for (let k = 0; k < 256; k++) img.data.set([lut[3 * k], lut[3 * k + 1], lut[3 * k + 2], 255], 4 * k);
  c.putImageData(img, 0, 0);
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

function panelImage(p, s) {
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

/** Draw a panel into `canvas` (w x h CSS px). Exports add a title and skip guides. */
function draw(p, canvas, w, h, dpr, exporting = false) {
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const c = canvas.getContext('2d');
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  c.fillStyle = '#ffffff';
  c.fillRect(0, 0, w, h);
  if (!p.data || !settings) return;

  const s = settings, X = meta.dims[p.x], Y = meta.dims[p.y], F = meta.dims[p.fixed];
  const g = geometry(p), sin = Math.sqrt(Math.max(0, 1 - g.cos * g.cos));
  const [u0, u1] = p.zoom?.u ?? viewRange(X, s.limit), [v0, v1] = p.zoom?.v ?? viewRange(Y, s.limit);
  const wx = (u, v) => g.lx * u + g.ly * g.cos * v, wy = (v) => g.ly * sin * v;
  const box = [[u0, v0], [u1, v0], [u1, v1], [u0, v1]];
  const xs = box.map(([u, v]) => wx(u, v)), ys = box.map(([, v]) => wy(v));
  const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys);
  const pad = { l: 62, r: 76, t: exporting ? 40 : 16, b: 52 };
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

  // Intensity image, clipped to the view parallelogram.
  const corners = box.map(([u, v]) => project(u, v));
  const { rows, cols } = p.data, ex = X.edges, ey = Y.edges;
  const dx = (ex[ex.length - 1] - ex[0]) / cols, dy = (ey[ey.length - 1] - ey[0]) / rows;
  c.save();
  c.beginPath();
  corners.forEach(([a, b], i) => (i ? c.lineTo(a, b) : c.moveTo(a, b)));
  c.closePath();
  c.fillStyle = MISSING;
  c.fill();
  c.clip();
  const [x0, y0] = project(ex[0], ey[0]);
  c.transform(sx * g.lx * dx, 0, sx * g.ly * g.cos * dy, -sy * g.ly * sin * dy, x0, y0);
  c.imageSmoothingEnabled = false;
  c.drawImage(panelImage(p, s), 0, 0);
  c.restore();

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
      for (const [color, offset] of [['rgba(255,255,255,.9)', 0], ['rgba(20,26,34,.55)', 4]]) {
        c.strokeStyle = color;
        c.lineDashOffset = offset;
        c.beginPath(); c.moveTo(a, b); c.lineTo(e, f); c.stroke();
      }
    }
    c.restore();
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
  for (const t of niceTicks(u0, u1)) {
    const [a, b] = project(t, v0);
    c.beginPath(); c.moveTo(a, b); c.lineTo(a, b + 4); c.stroke();
    c.fillText(fmt(t), a, b + 7);
  }
  c.textAlign = 'right';
  c.textBaseline = 'middle';
  for (const t of niceTicks(v0, v1)) {
    const [a, b] = project(u0, t);
    c.beginPath(); c.moveTo(a, b); c.lineTo(a - 4, b); c.stroke();
    c.fillText(fmt(t), a - 7, b);
  }
  c.fillStyle = INK;
  c.font = `600 12px ${SANS}`;
  c.textAlign = 'center';
  const [bx, by] = project((u0 + u1) / 2, v0);
  c.fillText(withUnits(X), bx, by + 34);
  const [ax0, ay0] = project(u0, v0), [ax1, ay1] = project(u0, v1);
  const len = Math.hypot(ax1 - ax0, ay1 - ay0), nx = (ay1 - ay0) / len, ny = -(ax1 - ax0) / len;
  c.save();
  c.translate((ax0 + ax1) / 2 + 46 * nx, (ay0 + ay1) / 2 + 46 * ny);
  c.rotate(Math.atan2(ay1 - ay0, ax1 - ax0));
  c.fillText(withUnits(Y), 0, 0);
  c.restore();

  if (exporting) {
    c.textAlign = 'left';
    c.textBaseline = 'alphabetic';
    c.font = `650 14px ${SANS}`;
    const title = `${F.label} = ${fmt(p.data.center)}`;
    c.fillText(title, 10, 24);
    const titleWidth = c.measureText(title).width;
    c.font = `11px ${MONO}`;
    c.fillStyle = INK2;
    const extras = `${symmetryNote(p.data)}${mask ? ` · mask ${pct((mask.edge + mask.outlier) / mask.measured)}` : ''}${p.data.removed ? ' · removed only' : ''}`;
    c.fillText(`slab ${fmt(p.data.slab[0])} to ${fmt(p.data.slab[1])}${extras}`, 20 + titleWidth, 24);
  }

  // Colorbar.
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

/** Redraw one panel, with the zoom box being dragged (as its u-v parallelogram). */
function drawPanel(p) {
  if (!settings) return;
  draw(p, p.canvas, p.canvas.clientWidth, p.canvas.clientHeight, devicePixelRatio || 1);
  const d = p.drag;
  if (!d?.moved || !p.project) return;
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
  for (const p of panels) drawPanel(p);
  paintColorbar();
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
  if (!pt) { p.hover.textContent = ''; return; }
  const [u, v] = pt, X = meta.dims[p.x], Y = meta.dims[p.y], ex = X.edges, ey = Y.edges;
  const { values, counts, rows, cols } = p.data;
  const col = Math.floor((u - ex[0]) / (ex[ex.length - 1] - ex[0]) * cols);
  const row = Math.floor((v - ey[0]) / (ey[ey.length - 1] - ey[0]) * rows);
  if (col < 0 || row < 0 || col >= cols || row >= rows) { p.hover.textContent = ''; return; }
  const i = row * cols + col, val = values[i], n = counts[i];
  p.hover.textContent = `${X.label} ${fmt(u)} ${Y.label} ${fmt(v)} → `
    + (Number.isFinite(val) ? `${fmtValue(val)} (${n} vox)` : 'no data');
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
  draw(p, out, p.canvas.clientWidth, p.canvas.clientHeight + 24, 3, true);
  const X = meta.dims[p.x], Y = meta.dims[p.y], F = meta.dims[p.fixed];
  out.toBlob((blob) => download(blob, `${stem()}_${X.label}${Y.label}_${F.label}=${p.data.center}.png`));
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
  showSymmetry();
  describe();
  panels.forEach(request);
  requestIso();
}

function showSymmetry() {
  const { ops } = symmetry, status = $('sym-status');
  status.className = 'note';
  $('sym-badge').textContent = ops.length > 1 ? `${symmetry.name} · ${ops.length}` : '1';
  if (ops.length === 1) {
    status.textContent = 'No symmetry averaging: each voxel is used as measured.';
  } else {
    let text = `${ops.length} operations; equivalent voxels are pooled with equal weight.`;
    const hkl = meta.dims.every((d) => d.basis);
    if (meta.lattice && hkl) {
      const change = metricChange(ops, meta.lattice);
      if (change > 0.02) {
        text = `${ops.length} operations, but they change the cell metric by up to ${(100 * change).toFixed(0)}%. `
          + 'They are not symmetries of this lattice; check the setting.';
        status.className = 'note warn';
      } else {
        text += ` Cell metric (${meta.lattice.source}) preserved to ${(100 * change).toFixed(2)}%.`;
      }
    } else if (!hkl) {
      text += ' Axes have no HKL basis, so operations act on the display axes directly.';
    }
    status.textContent = text;
  }
  $('sym-list').textContent = ops.map(formatOp).join('   ');
  $('sym-count').textContent = ops.length;
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
  engine('busy', 'building mask');
  $('mask-status').className = 'note';
  $('mask-status').textContent = radius || k ? 'Building mask…' : 'Clearing mask…';
  worker.postMessage({ type: 'mask', id: ++requestId, radius, k, maps: symmetry.maps, symmetry: symmetry.name });
}

function showMask(seconds) {
  const status = $('mask-status');
  status.className = 'note';
  if (!mask) {
    status.textContent = 'No mask: all measured voxels are used.';
    return;
  }
  const { measured, edge, outlier, radius, k, group } = mask;
  const parts = [];
  if (radius) parts.push(`${pct(edge / measured)} within ${radius} voxel${radius === 1 ? '' : 's'} of coverage edges`);
  if (k) parts.push(`${pct(outlier / measured)} above ${k}σ of their ${group} equivalents`);
  status.className = 'note ok';
  status.textContent = `Removed ${pct((edge + outlier) / measured)} of measured voxels: ${parts.join(', ')}${seconds ? ` (${seconds.toFixed(1)} s)` : ''}.`;
}

// ---- 3-D view ---------------------------------------------------------------------------

async function setup3D() {
  const section = document.createElement('section');
  section.className = 'card';
  section.innerHTML = `
    <div class="card-head">
      <span class="card-titlegroup"><span class="badge">3D</span><span class="panel-title">Isosurface and orthoslices</span>
        <span class="card-sub">drag to rotate · scroll to zoom · right-drag to pan</span></span>
      <span class="inline"><button type="button" class="btn btn-ghost btn-xs reset">Reset view</button>
        <button type="button" class="btn btn-ghost btn-xs save">Save PNG</button></span>
    </div>
    <canvas class="plot view3d" role="img" aria-label="3-D isosurface with the three slices"></canvas>
    <div class="controls3d">
      <div class="field">
        <div class="field-row"><span class="field-label">Isosurface level · log scale</span>
          <span class="readout-edit">level<input class="iso-level" type="number" step="any" style="width: 88px"></span></div>
        <input class="iso-slider" type="range" min="0" max="1000">
      </div>
      <div class="field">
        <div class="field-row"><span class="field-label">Surface opacity</span></div>
        <input class="iso-opacity" type="range" min="0.05" max="1" step="0.05" value="0.6">
      </div>
      <div class="field">
        <div class="field-row"><span class="field-label">Slice opacity</span></div>
        <input class="slice-opacity" type="range" min="0.05" max="1" step="0.05" value="1">
      </div>
      <div class="field"><span class="field-label">Grid</span>
        <div class="segmented iso-grid" data-value="100"><button type="button" data-value="64">Coarse</button><button type="button" data-value="100" class="on">Medium</button><button type="button" data-value="150">Fine</button></div></div>
      <label class="switch" style="padding-bottom: 6px"><input class="iso-slices" type="checkbox" checked><span class="switch-track"></span><span class="switch-label">Slices</span></label>
    </div>
    <div class="card-foot caption-foot"><span class="caption" aria-live="polite">Loading 3-D view…</span></div>`;
  $('view3d').append(section);
  const q = (sel) => section.querySelector(sel);
  iso = {
    levelInput: q('.iso-level'), slider: q('.iso-slider'), opacity: q('.iso-opacity'), sliceOpacity: q('.slice-opacity'), grid: q('.iso-grid'),
    slices: q('.iso-slices'), caption: q('.caption'), canvas: q('canvas'), userLevel: false, range: null,
  };
  paintAll();
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
  segmented(iso.grid, () => { iso.userLevel = false; requestIso(); });
  iso.slices.onchange = update3D;
  q('.reset').onclick = () => view3d.resetView();
  q('.save').onclick = () => view3d.snapshot((blob) => download(blob, `${stem()}_3d.png`));
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
  iso.wanted = { maxBins: Number(iso.grid.dataset.value), level, ops: symmetry.ops, symmetry: symmetry.name };
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
    const X = meta.dims[p.x], Y = meta.dims[p.y];
    return {
      fixed: p.fixed, x: p.x, y: p.y, center: p.data.center, image: panelImage(p, settings), key: p.imageKey,
      u: viewRange(X, settings.limit), v: viewRange(Y, settings.limit),
      ex: [X.edges[0], X.edges[X.edges.length - 1]], ey: [Y.edges[0], Y.edges[Y.edges.length - 1]],
    };
  });
  view3d.setSlices(slices, iso.slices.checked);
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
  try {
    const raw = localStorage.getItem('nxv-focus');
    const f = raw === null ? undefined : JSON.parse(raw);
    if (f === null || [0, 1, 2].includes(f)) focus = f;
  } catch { /* storage unavailable: keep the default */ }
  paintColorbar();
}

function persist() {
  const values = { cmap: $('cmap').value, scale: $('scale').dataset.value, angles: $('angles').checked, guides: $('guides').checked };
  try { localStorage.setItem('nxv-settings', JSON.stringify(values)); } catch { /* storage unavailable */ }
}

restore();
show('intro');
engine('', 'no file open');
$('open').onclick = $('open-intro').onclick = () => $('file').click();
$('file').onchange = () => { if ($('file').files[0]) openFile($('file').files[0]); $('file').value = ''; };
for (const id of ['cmap', 'vmin', 'vmax', 'soft', 'limit', 'angles', 'guides']) $(id).addEventListener('input', redraw);
for (const id of ['cmap', 'angles', 'guides']) $(id).addEventListener('change', persist);
$('cmap').addEventListener('input', paintColorbar);
segmented($('scale'), (value) => {
  if (value === 'log' && !(Number($('vmin').value) > 0)) $('vmin').value = sig(Number($('vmax').value) / 1000 || 1);
  persist();
  redraw();
});
$('auto').onclick = autoRange;
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
$('mask-removed').onchange = () => panels.forEach(request);
$('mask-download').onclick = () => worker.postMessage({ type: 'mask-download' });
const resized = new ResizeObserver(() => redraw());
resized.observe($('slices'));
new MutationObserver(() => { for (const c of $('slices').querySelectorAll('canvas')) resized.observe(c); }).observe($('slices'), { childList: true });

document.addEventListener('dragover', (e) => { e.preventDefault(); document.body.classList.add('dragging'); });
document.addEventListener('dragleave', (e) => { if (!e.relatedTarget) document.body.classList.remove('dragging'); });
document.addEventListener('drop', (e) => {
  e.preventDefault();
  document.body.classList.remove('dragging');
  const file = e.dataTransfer?.files?.[0];
  if (file) openFile(file);
});

const remote = new URLSearchParams(location.search).get('url');
if (remote) openURL(remote);
