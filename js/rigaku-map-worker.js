// Module worker for the parallel passes of reduceRigaku (rigaku-reduce.js), driven through a
// WorkerPool (rigaku-pool.js). It runs the same per-frame functions as the serial path; the
// coordinating worker merges the replies so that the result is bit-identical.
//   ping                          -> {}
//   sum-begin {nx, ny}            -> {}            start per-pixel sums (pass 1)
//   sum {buf, run, label}         -> {h}           decode, check and add one frame
//   sum-end                       -> {sum, counted, runSum}
//   peaks {bufs, hs, mask, bg}    -> {peaks}       one run's Bragg peaks (pass 2)
//   map-init {setup}              -> {}            pass 3 context (mapContext)
//   map {buf, h, run, grid, b}    -> {v, s, e, w, b}  one frame's contributions, grouped by voxel (pass 3);
//                                   b: the background counts of its unmasked pixels, or null
// Errors are returned as {error}.

import * as fmt from './rigaku-format.js';
import { frameRecords, indexFrame, makeSums, mapContext, mapFrame, runPeaks, sumFrame } from './rigaku-reduce.js';

let sums = null, img = null, ctx = null, frame = null, vox = null;

const handlers = {
  ping: () => ({}),
  'sum-begin': ({ nx, ny }) => {
    sums = makeSums(nx, ny);
    img = new Int32Array(nx * ny);
    return {};
  },
  sum: ({ buf, run, label }) => {
    const h = fmt.parseRodHeader(buf);
    sumFrame(fmt, buf, h, img, sums, run, label);
    return { h };
  },
  'sum-end': () => {
    const out = { sum: sums.sum, counted: sums.counted, runSum: sums.runSum };
    sums = null;
    return [out, [out.sum.buffer, out.counted.buffer, ...[...out.runSum.values()].map((a) => a.buffer)]];
  },
  peaks: ({ bufs, hs, mask, bg, nx, ny }) => ({ peaks: runPeaks(fmt, bufs, hs, mask, bg, nx, ny) }),
  'map-init': ({ setup }) => {
    ctx = mapContext(setup);
    frame = { x: new Float64Array(3 * ctx.np), c: new Int32Array(ctx.np), w: new Float64Array(ctx.np) };
    vox = new Int32Array(ctx.np);
    return {};
  },
  map: ({ buf, h, run, grid, b }) => {
    const A = mapFrame(ctx, fmt, buf, h, run, frame);
    const rec = frameRecords(vox, frame.c, frame.w, indexFrame(grid, frame.x, A, vox), ctx.nSub, grid.shape[0] * grid.shape[1] * grid.shape[2], b);
    return [rec, [rec.v.buffer, rec.s.buffer, rec.e.buffer, rec.w.buffer, ...(rec.b ? [rec.b.buffer] : [])]];
  },
};

self.onmessage = ({ data }) => {
  try {
    const r = handlers[data.type](data);
    const [out, transfer] = Array.isArray(r) ? r : [r, []];
    self.postMessage({ ...out, id: data.id }, transfer);
  } catch (err) {
    self.postMessage({ id: data.id, error: err?.message ?? String(err) });
  }
};
