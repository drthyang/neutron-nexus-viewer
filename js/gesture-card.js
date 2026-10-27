// The gesture card, after NEBULA3D's (web/src/components/GestureHint.tsx): the
// "?" beside the click modes opens a popover that shows how to zoom, move and
// reset a view and use the Box zoom and Cut tools, one row per action, the ways
// drawn as key caps and pointer glyphs (gestures.js, which keeps the rows true to
// app.js): the mouse, trackpad and keyboard ways with ⌘ or Ctrl, or the touch
// ones on a touch-only device, following that if it changes.
//
// A click or tap on the "?" (Enter or Space when it has focus) opens it and moves
// focus into it. Esc, a press outside it, the ×, or focus leaving it closes it;
// Esc, the × and a second press on the "?" hand focus back to the "?", while a
// press or Tab elsewhere leaves focus where it went. It does not open on hover:
// in a header of controls, a card that opens under a resting pointer covers the
// display bar. Until it has been opened once in this browser the "?" carries a
// small dot; nothing opens by itself.

import { gestureCard, isMacLike, isTouchOnly, platformString, TOUCH_QUERIES } from './gestures.js';

/** Where this browser remembers that the gesture card has been opened ("1"). */
export const SEEN_KEY = 'nxv-gesture-card-seen';

/** Apple platform: the modifier is ⌘. */
export const IS_MAC = isMacLike(platformString(navigator), navigator.userAgent);

const esc = (s) => String(s).replace(/[&<>"']/g, (ch) => `&#${ch.charCodeAt(0)};`);

/** A key cap. `name` is what a screen reader says for a symbol key ("Command" for ⌘). */
export function keyCap(k, name) {
  return name && name !== k
    ? `<kbd class="kbd"><span aria-hidden="true">${esc(k)}</span><span class="vh">${esc(name)}</span></kbd>`
    : `<kbd class="kbd">${esc(k)}</kbd>`;
}

/** The modifier's key cap: ⌘ on Apple platforms, Ctrl elsewhere. */
export const modKeyCap = (mac = IS_MAC) => (mac ? keyCap('⌘', 'Command') : keyCap('Ctrl'));

const svg = (body, size) => `<svg width="${size}" height="${size}" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${body}</svg>`;

// NEBULA3D's icons for zoom, move and reset; Box zoom and Cut have their header buttons' icons.
const ACTION_ICONS = {
  zoom: svg('<circle cx="7" cy="7" r="4.5"/><path d="M10.4 10.4 14 14M7 5v4M5 7h4"/>', 16),
  move: svg('<path d="M8 1.75v12.5M1.75 8h12.5M8 1.75 6.25 3.5M8 1.75 9.75 3.5M8 14.25 6.25 12.5M8 14.25l1.75-1.75M1.75 8 3.5 6.25M1.75 8 3.5 9.75M14.25 8 12.5 6.25M14.25 8 12.5 9.75"/>', 16),
  reset: svg('<path d="M2.75 8a5.25 5.25 0 1 0 1.6-3.8"/><path d="M2.5 2.5v3h3"/>', 16),
  box: svg('<path d="M1.75 4.25v-2.5h2.5M9.75 1.75h2.5v2.5M1.75 9.75v2.5h2.5M6.25 1.75h1.5M1.75 6.25v1.5M12.25 6.25v1M6.25 12.25h1"/><circle cx="10.75" cy="10.75" r="2.5"/><path d="M12.6 12.6 14.5 14.5"/>', 16),
  cut: svg('<path d="M3.5 12.5l9-9"/><circle cx="3" cy="13" r="1.6" fill="currentColor" stroke="none"/><circle cx="13" cy="3" r="1.6" fill="currentColor" stroke="none"/>', 16),
};

// The pointer glyphs, NEBULA3D's.
const GLYPHS = {
  scroll: svg('<rect x="4.5" y="1.75" width="7" height="12.5" rx="3.5"/><path d="M8 4.25v2.5"/>', 14), // a mouse with its wheel
  pinch: svg('<circle cx="4" cy="4" r="1.6"/><circle cx="12" cy="12" r="1.6"/><path d="m6.5 6.5 2 2M9.5 9.5l-1-1M8.5 6.5h-2v2M7.5 9.5h2v-2"/>', 14), // two fingertips drawn together
  drag: svg('<path d="M3 2.5v8l2.2-2 1.6 3.6 1.5-.7-1.6-3.5H9.7Z"/><path d="M11 11.5h3.5M12.8 9.8l1.7 1.7-1.7 1.7"/>', 14), // a pointer dragged
  'middle-drag': svg('<rect x="4.5" y="1.75" width="7" height="12.5" rx="3.5"/><path d="M8 3.75v3.5" stroke-width="2.6"/>', 14), // the middle button pressed
  'double-click': svg('<path d="M5 5.5v8l2.2-2 1.6 3.6 1.5-.7-1.6-3.5h2.4Z"/><path d="M3.2 3.2 2 2M6 2.4V.9M8.8 3.2 10 2"/>', 14), // a pointer with two clicks
  'two-finger': svg('<circle cx="4" cy="5" r="1.6"/><circle cx="4" cy="11" r="1.6"/><path d="M7.5 8h6.5M12.3 6.3 14 8l-1.7 1.7"/>', 14), // two fingertips moved together
  'double-tap': svg('<circle cx="8" cy="8" r="1.8"/><circle cx="8" cy="8" r="4" opacity="0.7"/><circle cx="8" cy="8" r="6.3" opacity="0.4"/>', 14), // a fingertip tapped twice
};

function tokenHTML(tk) {
  switch (tk.t) {
    case 'key': return keyCap(tk.k, tk.name);
    case 'button': return `<span class="gh-btn">${esc(tk.label)}</span>`;
    case 'glyph': return `<span class="gh-glyph">${GLYPHS[tk.g]}<span>${esc(tk.text)}</span></span>`;
    default: return `<span class="gh-text">${esc(tk.text)}</span>`;
  }
}

// The ways of a row, with a comma between them for a screen reader.
const rowHTML = (r) => `<li class="gh-row" data-action="${esc(r.id)}"><span class="gh-icon">${ACTION_ICONS[r.id]}</span><span class="gh-name">${esc(r.name)}</span><span class="gh-ways">${
  r.ways.map((w, i) => `${i ? '<span class="vh">, </span>' : ''}<span class="gh-way">${w.map(tokenHTML).join('')}</span>`).join('')}</span></li>`;

/**
 * Wire the gesture card in `wrap` (index.html's #gesture-hint: the "?" button,
 * .gh-trigger, and the card, .gh-pop > .gh-card, whose rows are drawn here).
 */
export function setupGestureCard(wrap) {
  const trigger = wrap.querySelector('.gh-trigger'), pop = wrap.querySelector('.gh-pop'), card = wrap.querySelector('.gh-card');
  const sub = card.querySelector('.gh-sub'), rows = card.querySelector('.gh-rows'), note = card.querySelector('.gh-note');
  const media = typeof matchMedia === 'function' ? TOUCH_QUERIES.map((q) => matchMedia(q)) : [];
  const touchOnly = () => media.length > 0 && isTouchOnly((q) => matchMedia(q).matches);
  const isOpen = () => !pop.hidden;
  // A press inside the card or on the "?" that has not lifted yet: Safari moves
  // focus nowhere on a button press, which must not count as focus leaving.
  let pressing = false;

  let seen = false;
  try { seen = localStorage.getItem(SEEN_KEY) === '1'; } catch { /* storage unavailable: the dot shows until the card is opened this visit */ }
  trigger.classList.toggle('gh-new', !seen);

  function render() {
    const c = gestureCard({ mac: IS_MAC, touch: touchOnly() });
    sub.textContent = c.sub;
    rows.innerHTML = c.rows.map(rowHTML).join('');
    note.textContent = c.note ?? '';
    note.hidden = !c.note;
  }

  // Keep the card inside the window: shift it left when the "?" sits near the
  // right edge, narrow it on a narrow screen, and let it scroll in a short window.
  function place() {
    const vw = document.documentElement.clientWidth;
    card.style.maxWidth = `${vw - 16}px`;
    card.style.maxHeight = '';
    const left = wrap.getBoundingClientRect().left, width = card.offsetWidth;
    pop.style.left = `${Math.min(Math.max(8, left), vw - width - 8) - left}px`;
    card.style.maxHeight = `${Math.max(160, innerHeight - card.getBoundingClientRect().top - 8)}px`;
  }

  function open() {
    render();
    pop.hidden = false;
    trigger.setAttribute('aria-expanded', 'true');
    trigger.setAttribute('aria-controls', card.id);
    place();
    card.focus();
    // Opened once: the dot goes, for good.
    if (!seen) {
      seen = true;
      trigger.classList.remove('gh-new');
      try { localStorage.setItem(SEEN_KEY, '1'); } catch { /* storage unavailable: seen for this visit only */ }
    }
    addEventListener('keydown', onKey, true);
    document.addEventListener('pointerdown', onDown, true);
    addEventListener('resize', place);
  }

  function close(refocus) {
    if (!isOpen()) return;
    pop.hidden = true;
    trigger.setAttribute('aria-expanded', 'false');
    trigger.removeAttribute('aria-controls');
    removeEventListener('keydown', onKey, true);
    document.removeEventListener('pointerdown', onDown, true);
    removeEventListener('resize', place);
    if (refocus) trigger.focus();
  }

  // Esc closes it, before the page's own Esc (the other popovers, an armed Box zoom, the single layout).
  function onKey(e) {
    if (e.key !== 'Escape') return;
    e.preventDefault();
    e.stopPropagation();
    close(true);
  }
  // A press outside closes it.
  function onDown(e) {
    if (!wrap.contains(e.target)) close(false);
  }

  trigger.addEventListener('click', () => (isOpen() ? close(true) : open()));
  card.querySelector('.gh-close').addEventListener('click', () => close(true));
  wrap.addEventListener('pointerdown', () => { pressing = true; });
  document.addEventListener('pointerup', () => { pressing = false; }, true);
  document.addEventListener('pointercancel', () => { pressing = false; }, true);
  // Focus leaving it (Tab past its end, or to another control) closes it.
  wrap.addEventListener('focusout', (e) => {
    if (isOpen() && !pressing && !wrap.contains(e.relatedTarget)) close(false);
  });
  // A keyboard and trackpad attached to a tablet, or taken away: the other rows.
  for (const m of media) m.addEventListener?.('change', () => { if (isOpen()) { render(); place(); } });
  return { open, close, isOpen };
}
