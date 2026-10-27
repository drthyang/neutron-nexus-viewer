// Zoom and pan of a view's window, for the slices and the 1-D plots (pure
// functions). The zoom rate and the wheel rules are NEBULA3D's (its slice
// viewers' viewport.ts and ViewFrame.tsx), so the two apps zoom alike.
//
// A slice's window is { u: [u0, u1], v: [v0, v1] } in its axes' units; a plot's
// x range is [x0, x1]. Every zoom and pan keeps the window's centre on the data,
// as a map keeps its centre on the world: however far it is panned, the data
// stays on screen.

/**
 * The zoom factor of one wheel event (< 1 zooms in): its vertical delta read
 * in pixels whatever the device reports (a line as 40 px, a page as the view's
 * height), e^(Δ/100 px), at most 2× either way. One mouse notch (100 px in
 * Chrome, 3 lines in Firefox) is then a factor 2, and a trackpad's or a
 * pinch's small steps zoom smoothly.
 */
export function wheelZoomFactor(deltaY, deltaMode, viewHeight) {
  const px = deltaMode === 1 ? deltaY * 40 : deltaMode === 2 ? deltaY * viewHeight : deltaY;
  return Math.min(2, Math.max(0.5, Math.exp(px * 0.01)));
}

/** Move a range so its centre lies in `extent` = [lo, hi]; a range whose centre does is returned as it is. */
export function keepCentreIn(range, extent) {
  const c = (range[0] + range[1]) / 2, m = Math.min(extent[1], Math.max(extent[0], c));
  return m === c ? range : [range[0] + m - c, range[1] + m - c];
}

const scaleAbout = (range, at, k) => [at + (range[0] - at) * k, at + (range[1] - at) * k];
// Zooming out to (nearly) the full view's size is the full view.
const reaches = (size, full) => size >= full * (1 - 1e-9);

/**
 * Zoom a slice window by `k` (< 1 zooms in) about the point `at` = [u, v],
 * which keeps its place on screen. Both axes scale alike, so the window keeps
 * its shape. `limits`:
 *   least   [du, dv], the smallest size on each axis (two bins): k stops there
 *   full    the full view's window: zooming out to its size on either axis
 *           returns null, the full view
 *   extent  { u: [lo, hi], v: [lo, hi] }, the data, where the centre stays
 */
export function zoomWindow(win, at, k, { least, full, extent }) {
  const su = win.u[1] - win.u[0], sv = win.v[1] - win.v[0];
  // Not below `least` on either axis; a window already smaller only zooms out.
  k = Math.max(k, Math.min(1, Math.max(least[0] / su, least[1] / sv)));
  if (k > 1 && (reaches(su * k, full.u[1] - full.u[0]) || reaches(sv * k, full.v[1] - full.v[0]))) return null;
  return { u: keepCentreIn(scaleAbout(win.u, at[0], k), extent.u), v: keepCentreIn(scaleAbout(win.v, at[1], k), extent.v) };
}

/** Shift a slice window by (du, dv), keeping its centre on the data (`extent`, as for zoomWindow()). */
export function panWindow(win, du, dv, extent) {
  return { u: keepCentreIn([win.u[0] + du, win.u[1] + du], extent.u), v: keepCentreIn([win.v[0] + dv, win.v[1] + dv], extent.v) };
}

/**
 * Zoom a plot's range by `k` (< 1 zooms in) about `at`: at least `least`
 * across, its centre in `extent` (default `full`), and null, the full range,
 * once zooming out reaches `full`'s size.
 */
export function zoomRange(range, at, k, { least, full, extent = full }) {
  const size = range[1] - range[0];
  k = Math.max(k, Math.min(1, least / size));
  if (k > 1 && reaches(size * k, full[1] - full[0])) return null;
  return keepCentreIn(scaleAbout(range, at, k), extent);
}

/**
 * A plot's range under two fingers, from `range` at the start of the gesture:
 * scaled by `k` (the fingers' spread at the start over now; < 1 zooms in) and
 * placed so `at`, the value under their midpoint at the start, is under it now,
 * a fraction `f` along the range. As in a slice's pinch (and NEBULA3D's), it
 * zooms and pans in one, between `least` across and `full`'s size, with its
 * centre in `extent` (default `full`). It does not snap to the full range
 * mid-gesture, so a two-finger pan at the full size does not flicker.
 */
export function pinchRange(range, at, f, k, { least, full, extent = full }) {
  const size0 = range[1] - range[0];
  k = Math.max(k, Math.min(1, least / size0));
  const size = Math.min(size0 * k, Math.max(size0, full[1] - full[0]));
  const lo = at - f * size;
  return keepCentreIn([lo, lo + size], extent);
}
