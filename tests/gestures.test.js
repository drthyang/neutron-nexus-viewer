import assert from 'node:assert/strict';
import test from 'node:test';

import {
  DBLCLICK_AFTER_TAP_MS, DOUBLE_TAP_MS, DOUBLE_TAP_PX, TAP_SLOP,
  gestureCard, isDoubleTap, isMacLike, isTouchOnly, modKey, platformString, tapTracker, wayText,
} from '../js/gestures.js';

const ways = (card, id) => card.rows.find((r) => r.id === id).ways.map(wayText);

test('⌘ on Apple platforms, iPadOS (which reports MacIntel) included; Ctrl elsewhere', () => {
  for (const p of ['MacIntel', 'macOS', 'iPhone', 'iPad', 'iPod touch']) assert.ok(isMacLike(p), p);
  for (const p of ['Win32', 'Windows', 'Linux x86_64', 'Linux armv8l', 'Android', 'CrOS']) assert.ok(!isMacLike(p), p);
  // No platform: the user agent decides.
  assert.ok(isMacLike('', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'));
  assert.ok(!isMacLike('', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'));
  assert.deepEqual(modKey(true), { t: 'key', k: '⌘', name: 'Command' });
  assert.deepEqual(modKey(false), { t: 'key', k: 'Ctrl' });
});

test('the platform comes from Client Hints where a browser has them, else navigator.platform', () => {
  assert.equal(platformString({ userAgentData: { platform: 'macOS' }, platform: 'MacIntel' }), 'macOS');
  assert.equal(platformString({ userAgentData: { platform: '' }, platform: 'Win32' }), 'Win32');
  assert.equal(platformString({ platform: 'iPad' }), 'iPad');
  assert.equal(platformString(undefined), '');
});

test('touch-only means a coarse pointer and no hover', () => {
  const media = (on) => (q) => on.includes(q);
  assert.ok(isTouchOnly(media(['(pointer: coarse)', '(hover: none)'])));
  assert.ok(!isTouchOnly(media(['(pointer: coarse)']))); // a touch laptop: it hovers
  assert.ok(!isTouchOnly(media(['(hover: none)'])));
  assert.ok(!isTouchOnly(media([])));
});

test('the desktop card: the mouse, trackpad and keyboard ways, with the platform\'s modifier', () => {
  const mac = gestureCard({ mac: true, touch: false }), pc = gestureCard({ mac: false, touch: false });
  assert.deepEqual(mac.rows.map((r) => r.name), ['Zoom', 'Move', 'Reset', 'Box zoom', 'Cut']);
  assert.deepEqual(ways(mac, 'zoom'), ['Command + scroll', 'pinch', '+ minus']);
  assert.deepEqual(ways(pc, 'zoom'), ['Ctrl + scroll', 'pinch', '+ minus']);
  assert.deepEqual(ways(mac, 'move'), ['Shift + drag', 'middle-drag', 'Left Up Down Right']);
  assert.deepEqual(ways(mac, 'reset'), ['double-click', '0', 'Reset zoom']);
  assert.deepEqual(ways(mac, 'box'), ['Box zoom then drag a box', 'Esc cancels']);
  assert.deepEqual(ways(mac, 'cut'), ['Cut then drag across a slice', 'Shift mid-drag: lattice direction']);
  assert.equal(mac.sub, 'in any mode, on the slices and plots');
  assert.equal(mac.note, 'Keys act on the view you clicked last.');
  // Only the modifier differs.
  assert.deepEqual(pc.rows.slice(1), mac.rows.slice(1));
});

test('the touch card: pinch, two-finger drag and double-tap, and no keys', () => {
  const card = gestureCard({ mac: true, touch: true });
  assert.deepEqual(card.rows.map((r) => r.name), ['Zoom', 'Move', 'Reset', 'Box zoom', 'Cut']);
  assert.deepEqual(ways(card, 'zoom'), ['pinch']);
  assert.deepEqual(ways(card, 'move'), ['two-finger drag']);
  assert.deepEqual(ways(card, 'reset'), ['double-tap', 'Reset zoom']);
  assert.deepEqual(ways(card, 'box'), ['Box zoom then drag a box', 'tap it again to cancel']);
  assert.deepEqual(ways(card, 'cut'), ['Cut then drag across a slice']);
  assert.ok(card.rows.every((r) => r.ways.flat().every((tk) => tk.t !== 'key')));
  assert.equal(card.note, undefined);
  assert.deepEqual(gestureCard({ mac: false, touch: true }), card);
});

test('a double-tap: two taps within 320 ms and 30 px', () => {
  assert.equal(DOUBLE_TAP_MS, 320);
  assert.equal(DOUBLE_TAP_PX, 30);
  const a = { t: 1000, x: 100, y: 100 };
  assert.ok(isDoubleTap(a, { t: 1300, x: 120, y: 110 }));
  assert.ok(!isDoubleTap(a, { t: 1320, x: 100, y: 100 })); // too late
  assert.ok(!isDoubleTap(a, { t: 1100, x: 130, y: 100 })); // too far
  assert.ok(!isDoubleTap(null, a));
});

test('the tap counter: taps, not drags or pinches, and one reset per double-tap', () => {
  const taps = tapTracker();
  taps.down(1, 100, 100);
  assert.equal(taps.up(1, 101, 100, 0), false); // the first tap
  taps.down(2, 110, 105);
  assert.equal(taps.up(2, 110, 105, 200), true); // the second: a double-tap
  // A third tap right after starts over rather than resetting again.
  taps.down(3, 110, 105);
  assert.equal(taps.up(3, 110, 105, 300), false);
  // The dblclick a browser may send after it is the same gesture.
  assert.ok(taps.justDoubleTapped(200 + DBLCLICK_AFTER_TAP_MS - 1));
  assert.ok(!taps.justDoubleTapped(200 + DBLCLICK_AFTER_TAP_MS));

  // A touch that strays past the slop is a drag, even when it comes back.
  const drag = tapTracker();
  drag.down(1, 50, 50);
  drag.up(1, 50, 50, 0);
  drag.down(2, 50, 50);
  drag.move(2, 50 + TAP_SLOP + 1, 50);
  drag.move(2, 50, 50);
  assert.equal(drag.up(2, 50, 50, 100), false);

  // The fingers of a pinch are not taps, and a tap after one is a first tap.
  const pinch = tapTracker();
  pinch.down(1, 50, 50);
  pinch.up(1, 50, 50, 0);
  pinch.down(2, 50, 50);
  pinch.down(3, 150, 150);
  assert.equal(pinch.up(2, 50, 50, 100), false);
  assert.equal(pinch.up(3, 150, 150, 120), false);
  pinch.down(4, 50, 50);
  assert.equal(pinch.up(4, 50, 50, 200), false);
  pinch.down(5, 50, 50);
  assert.equal(pinch.up(5, 50, 50, 300), true);

  // A cancelled touch forgets the tap before it; an unknown pointer is ignored.
  const cancelled = tapTracker();
  cancelled.down(1, 0, 0);
  cancelled.up(1, 0, 0, 0);
  cancelled.down(2, 0, 0);
  cancelled.cancel(2);
  cancelled.down(3, 0, 0);
  assert.equal(cancelled.up(3, 0, 0, 100), false);
  assert.equal(cancelled.up(9, 0, 0, 150), false);
});
