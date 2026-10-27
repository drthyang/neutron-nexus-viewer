import assert from 'node:assert/strict';
import test from 'node:test';

import { keepCentreIn, panWindow, pinchRange, wheelZoomFactor, zoomRange, zoomWindow } from '../js/viewzoom.js';

const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) < eps, `${a} ≠ ${b}`);
const closeAll = (a, b) => a.forEach((x, i) => close(x, b[i]));

test('the wheel zoom rate is NEBULA3D\'s: one mouse notch is 2×, small steps are smooth', () => {
  // Chrome's notch is 100 px, Firefox's 3 lines (40 px each): both 2× at most.
  close(wheelZoomFactor(100, 0, 500), 2);
  close(wheelZoomFactor(-100, 0, 500), 0.5);
  close(wheelZoomFactor(3, 1, 500), 2);
  close(wheelZoomFactor(-3, 1, 500), 0.5);
  // A trackpad or pinch step of a few pixels.
  close(wheelZoomFactor(-4, 0, 500), Math.exp(-0.04));
  // A page is the view's height, capped at 2×.
  close(wheelZoomFactor(1, 2, 30), Math.exp(0.3));
  close(wheelZoomFactor(1, 2, 500), 2);
});

test('a range keeps its centre on the data', () => {
  assert.deepEqual(keepCentreIn([1, 3], [0, 10]), [1, 3]);
  assert.deepEqual(keepCentreIn([-4, -2], [0, 10]), [-1, 1]);
  assert.deepEqual(keepCentreIn([11, 13], [0, 10]), [9, 11]);
});

const limits = { least: [0.2, 0.2], full: { u: [-5, 5], v: [-4, 4] }, extent: { u: [-5, 5], v: [-4, 4] } };

test('zooming a slice window keeps the point under the pointer and the window\'s shape', () => {
  const win = { u: [-5, 5], v: [-4, 4] };
  const z = zoomWindow(win, [2, 1], 0.5, limits);
  closeAll(z.u, [-1.5, 3.5]);
  closeAll(z.v, [-1.5, 2.5]);
  // The anchor sits at the same fraction of the window on both axes.
  close((2 - z.u[0]) / (z.u[1] - z.u[0]), (2 - win.u[0]) / 10);
  close((1 - z.v[0]) / (z.v[1] - z.v[0]), (1 - win.v[0]) / 8);
  // Zooming back out by the inverse factor returns to the full view.
  assert.equal(zoomWindow(z, [2, 1], 2, limits), null);
});

test('zooming stops at two bins and returns to the full view when zoomed out to its size', () => {
  const small = zoomWindow({ u: [0, 1], v: [0, 0.8] }, [0.5, 0.4], 0.01, limits);
  close(small.v[1] - small.v[0], 0.2);
  close(small.u[1] - small.u[0], 0.25);
  // At the limit already: the window stays.
  const again = zoomWindow(small, [0.5, 0.4], 0.5, limits);
  closeAll(again.u, small.u);
  // Zooming out stops being a zoom once either axis spans the full view.
  assert.equal(zoomWindow({ u: [-2, 2], v: [-3.9, 3.9] }, [0, 0], 1.1, limits), null);
  const out = zoomWindow({ u: [-2, 2], v: [-1, 1] }, [0, 0], 1.5, limits);
  closeAll(out.u, [-3, 3]);
  // From the full view, zooming out does nothing.
  assert.equal(zoomWindow(limits.full, [0, 0], 2, limits), null);
});

test('zooming and panning keep the window\'s centre on the data', () => {
  // Zooming in about a point off the data moves the centre back onto its edge.
  const z = zoomWindow({ u: [-5, 5], v: [-4, 4] }, [-20, 0], 0.5, { ...limits, full: { u: [-50, 50], v: [-40, 40] } });
  close((z.u[0] + z.u[1]) / 2, -5);
  close(z.u[1] - z.u[0], 5);
  const p = panWindow({ u: [0, 2], v: [0, 2] }, 100, -100, limits.extent);
  closeAll(p.u, [4, 6]);
  closeAll(p.v, [-5, -3]);
  const q = panWindow({ u: [0, 2], v: [0, 2] }, 1, -0.5, limits.extent);
  closeAll(q.u, [1, 3]);
  closeAll(q.v, [-0.5, 1.5]);
});

test('a plot\'s range zooms about the pointer, at least two bins across', () => {
  const opts = { least: 0.1, full: [0, 10] };
  closeAll(zoomRange([0, 10], 4, 0.5, opts), [2, 7]);
  closeAll(zoomRange([4, 6], 5, 0.001, opts), [4.95, 5.05]);
  assert.equal(zoomRange([2, 7], 4, 2, opts), null);
  // Its centre stays on the data.
  closeAll(zoomRange([0, 10], -10, 0.5, opts), [-2.5, 2.5]);
});

test('two fingers on a plot zoom about their midpoint and pan with it', () => {
  const opts = { least: 0.1, full: [0, 10] };
  // Spread twice as far, midpoint still: the value under it (4, 40% along) stays there.
  closeAll(pinchRange([0, 10], 4, 0.4, 0.5, opts), [2, 7]);
  // The same spread, the midpoint moved from 40% to 60% along: a pan by a fifth of the range.
  closeAll(pinchRange([2, 7], 4, 0.6, 1, opts), [1, 6]);
  // A pinch and a pan at once: the value 4 lands 25% along a range half as wide.
  closeAll(pinchRange([2, 7], 4, 0.25, 0.5, opts), [3.375, 5.875]);
  // Not below two bins across, nor out past the full range's size.
  closeAll(pinchRange([4, 6], 5, 0.5, 0.001, opts), [4.95, 5.05]);
  closeAll(pinchRange([2, 7], 4.5, 0.5, 8, opts), [-0.5, 9.5]); // the full size, 4.5 still mid-way
  // No snap to the full range mid-gesture: at the full size, a two-finger pan pans.
  closeAll(pinchRange([0, 10], 5, 0.6, 1.02, opts), [-1, 9]);
  // Its centre stays on the data.
  closeAll(pinchRange([0, 10], 5, -2, 1, opts), [5, 15]);
  closeAll(pinchRange([2, 4], 3, -10, 1, opts), [9, 11]);
});
