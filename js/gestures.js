// What the gesture card says (gesture-card.js draws it), after NEBULA3D's gesture
// card (web/src/components/gestures.ts) so the two apps read the same: one row per
// action, each with the ways to do it, drawn as key caps and pointer glyphs.
// Every way here is one app.js handles; keep the two in step:
//   zoom   Ctrl/⌘ + wheel or a pinch (wireZoom: the wheel and Safari's
//          GestureEvents; movePinch and moveCurvePinch: two fingers on a
//          slice or a plot);
//          + / − on a focused view (wireZoom's keys)
//   move   Shift + drag or a middle-button drag (panStart); the arrow keys on a
//          focused view; two fingers dragged together on a slice or a plot
//   reset  double-click (and a double-tap on a touch screen, tapTracker below),
//          0 on a focused view, Reset zoom in the view's header
//   box    the Box zoom button, then one drag; Esc or a second press puts it away
//   cut    the Cut button, then a drag across a slice; Shift pressed during the
//          drag keeps the cut along a lattice direction (moveCutDrag), while
//          Shift held before it pans
// The slices and the I(Q) and cut plots take all of these (the 3-D view has its
// own camera), the touch gestures too.
// Touch-only devices get the touch gestures; everything else the mouse,
// trackpad and keyboard ones, with ⌘ on Apple platforms and Ctrl elsewhere.
// Pure functions: no DOM here.

/** Apple platforms, whose modifier is ⌘ (iPadOS reports "MacIntel"). */
export function isMacLike(platform, userAgent = '') {
  return /mac|iphone|ipad|ipod/i.test(platform || userAgent);
}

/** The platform a browser reports: User-Agent Client Hints where there are any, else navigator.platform. */
export function platformString(nav) {
  return nav?.userAgentData?.platform || nav?.platform || '';
}

/** A touch screen with no hover: no mouse or trackpad to point with. */
export function isTouchOnly(matches) {
  return matches('(pointer: coarse)') && matches('(hover: none)');
}

/** The media queries isTouchOnly() reads, to follow their changes (a tablet's keyboard and trackpad attached). */
export const TOUCH_QUERIES = ['(pointer: coarse)', '(hover: none)'];

const key = (k, name) => ({ t: 'key', k, ...(name ? { name } : {}) }); // `name` is read out when `k` is a symbol
const button = (label) => ({ t: 'button', label }); // an on-screen button, drawn like one
const glyph = (g, text) => ({ t: 'glyph', g, text }); // a pointer gesture: its glyph and words
const text = (t) => ({ t: 'text', text: t });

/** The modifier key cap: ⌘ on Apple platforms, Ctrl elsewhere. */
export const modKey = (mac) => (mac ? key('⌘', 'Command') : key('Ctrl'));

/**
 * The card for this device: { sub, rows: [{ id, name, ways }], note? }, where a
 * way is a list of tokens ({ t: 'key' | 'button' | 'glyph' | 'text', … }).
 */
export function gestureCard({ mac, touch }) {
  const box = [button('Box zoom'), text('then'), glyph('drag', 'drag a box')];
  const cut = [button('Cut'), text('then'), glyph('drag', 'drag across a slice')];
  if (touch) {
    return {
      sub: 'in any mode, on the slices and plots',
      rows: [
        { id: 'zoom', name: 'Zoom', ways: [[glyph('pinch', 'pinch')]] },
        { id: 'move', name: 'Move', ways: [[glyph('two-finger', 'two-finger drag')]] },
        { id: 'reset', name: 'Reset', ways: [[glyph('double-tap', 'double-tap')], [button('Reset zoom')]] },
        { id: 'box', name: 'Box zoom', ways: [box, [text('tap it again to cancel')]] },
        { id: 'cut', name: 'Cut', ways: [cut] },
      ],
    };
  }
  return {
    sub: 'in any mode, on the slices and plots',
    rows: [
      { id: 'zoom', name: 'Zoom', ways: [[modKey(mac), text('+'), glyph('scroll', 'scroll')], [glyph('pinch', 'pinch')], [key('+'), key('−', 'minus')]] },
      {
        id: 'move',
        name: 'Move',
        ways: [
          [key('⇧ Shift', 'Shift'), text('+'), glyph('drag', 'drag')],
          [glyph('middle-drag', 'middle-drag')],
          [key('←', 'Left'), key('↑', 'Up'), key('↓', 'Down'), key('→', 'Right')],
        ],
      },
      { id: 'reset', name: 'Reset', ways: [[glyph('double-click', 'double-click')], [key('0')], [button('Reset zoom')]] },
      { id: 'box', name: 'Box zoom', ways: [box, [key('Esc'), text('cancels')]] },
      { id: 'cut', name: 'Cut', ways: [cut, [key('⇧ Shift', 'Shift'), text('mid-drag: lattice direction')]] },
    ],
    note: 'Keys act on the focused view (click it or Tab to it).',
  };
}

/** A way in words, for a screen reader or a test ("⌘ + scroll" → "Command + scroll"). */
export function wayText(way) {
  return way.map((tk) => (tk.t === 'key' ? tk.name ?? tk.k : tk.t === 'button' ? tk.label : tk.text)).join(' ');
}

// ---- Double-tap ------------------------------------------------------------------------
//
// Browsers do not reliably send dblclick for a double-tap, so the views count
// taps themselves, as NEBULA3D's do: two taps this close in time and place are a
// double-tap, which resets the view as a double-click does.

/** The most time between the two taps of a double-tap, in ms. */
export const DOUBLE_TAP_MS = 320;
/** The farthest apart the two taps of a double-tap may land, in CSS px. */
export const DOUBLE_TAP_PX = 30;
/** How far a touch may move and still be a tap, in CSS px (a slice's drag starts past 5 px). */
export const TAP_SLOP = 5;
/** A dblclick this soon after a double-tap is the same gesture, in ms: not a second reset. */
export const DBLCLICK_AFTER_TAP_MS = 600;

/** Whether tap `b` ({ t, x, y }) follows tap `a` closely enough to make a double-tap. */
export function isDoubleTap(a, b) {
  return !!a && b.t - a.t < DOUBLE_TAP_MS && Math.hypot(b.x - a.x, b.y - a.y) < DOUBLE_TAP_PX;
}

/**
 * Counts the taps of the touches on one view. Feed it the touch pointers'
 * down, move, up and cancel (id, CSS px, ms); up() says whether that touch
 * ended a double-tap. A tap is a touch that lifts within TAP_SLOP of where it
 * landed, never having strayed further, with no other finger down meanwhile
 * (a pinch is not a tap).
 */
export function tapTracker() {
  const down = new Map();
  let multi = false, last = null, resetAt = -Infinity;
  const lift = (id) => {
    down.delete(id);
    if (!down.size) multi = false;
  };
  return {
    down(id, x, y) {
      down.set(id, { x, y, moved: false });
      if (down.size > 1) { multi = true; last = null; }
    },
    move(id, x, y) {
      const d = down.get(id);
      if (d && Math.hypot(x - d.x, y - d.y) > TAP_SLOP) d.moved = true;
    },
    up(id, x, y, t) {
      const d = down.get(id);
      if (!d) return false;
      const tap = !multi && !d.moved && Math.hypot(x - d.x, y - d.y) <= TAP_SLOP;
      lift(id);
      if (!tap) { last = null; return false; }
      const here = { t, x, y };
      if (isDoubleTap(last, here)) { last = null; resetAt = t; return true; }
      last = here;
      return false;
    },
    cancel(id) { lift(id); last = null; },
    /** Whether a dblclick at time t belongs to the double-tap just counted. */
    justDoubleTapped: (t) => t - resetAt < DBLCLICK_AFTER_TAP_MS,
  };
}
