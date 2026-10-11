// Slices as PNG images for a model to look at: drawn in their true geometry (oblique axes at their
// real angle, equal lengths in reciprocal space equal on screen, as the viewer draws them), with the
// viewer's color scales, and no smoothing: every screen pixel shows the bin it falls in.

import { deflateSync } from 'node:zlib';

import { COLORMAPS } from '../js/colormaps.js';

const LUTS = new Map();
function lut(name) {
  if (!LUTS.has(name)) {
    const hex = COLORMAPS[name];
    if (!hex) throw new Error(`Unknown colormap "${name}".`);
    LUTS.set(name, Uint8Array.from({ length: 768 }, (_, i) => parseInt(hex.slice(2 * i, 2 * i + 2), 16)));
  }
  return LUTS.get(name);
}

export const COLORMAP_NAMES = Object.keys(COLORMAPS);

/** The viewer's mapping of a value to [0, 1] for scale 'asinh' (with softening), 'linear' or 'log'. */
export function colorScale({ scale, min, max, soft }) {
  const f = scale === 'asinh' ? (v) => Math.asinh(v / soft) : scale === 'log' ? (v) => Math.log10(Math.max(v, min)) : (v) => v;
  const lo = f(min), span = f(max) - lo;
  return (v) => Math.min(1, Math.max(0, (f(v) - lo) / (span || 1)));
}

/**
 * Render slice values (row-major, `rows` × `cols`; column = x bin, row = y bin) as RGB. `xEdges` and
 * `yEdges` are the bin edges, `geometry` the plane's { lx, ly, cos, equal } (planeGeometry()), so y
 * points up and meets x at its true angle. Bins without data are grey. Returns { png, width, height }.
 */
export function renderSlice({ values, rows, cols, xEdges, yEdges, geometry, color, colormap, maxSize = 640 }) {
  const { lx, ly, cos, equal } = geometry, sin = Math.sqrt(Math.max(0, 1 - cos * cos));
  const [x0, x1] = [xEdges[0], xEdges[cols]], [y0, y1] = [yEdges[0], yEdges[rows]];
  // Screen coordinates of a display point (u, v): X = u lx + v ly cos, Y = v ly sin (up).
  let ax = lx, bx = ly * cos, by = ly * sin;
  if (!equal) {
    // Axes in different units: fill a square instead.
    ax = 1 / (x1 - x0); bx = 0; by = 1 / (y1 - y0);
  }
  const corners = [[x0, y0], [x1, y0], [x0, y1], [x1, y1]].map(([u, v]) => [u * ax + v * bx, v * by]);
  const minX = Math.min(...corners.map((c) => c[0])), maxX = Math.max(...corners.map((c) => c[0]));
  const minY = Math.min(...corners.map((c) => c[1])), maxY = Math.max(...corners.map((c) => c[1]));
  // At most maxSize pixels on the long side, and no more than 4 pixels per bin.
  const binSize = Math.max(Math.abs(ax) * (x1 - x0) / cols, by * (y1 - y0) / rows);
  const scale = Math.min(maxSize / Math.max(maxX - minX, maxY - minY), 4 / binSize);
  const width = Math.max(1, Math.round((maxX - minX) * scale)), height = Math.max(1, Math.round((maxY - minY) * scale));
  const table = lut(colormap), map = colorScale(color);
  const xFind = finder(xEdges), yFind = finder(yEdges);
  const rgb = new Uint8Array(width * height * 3).fill(255);
  for (let py = 0; py < height; py++) {
    const Y = maxY - (py + 0.5) / scale, v = Y / by;
    const row = yFind(v);
    for (let px = 0; px < width; px++) {
      const X = minX + (px + 0.5) / scale, u = (X - v * bx) / ax;
      const col = xFind(u);
      if (row < 0 || col < 0) continue; // outside the plane: white
      const value = values[row * cols + col], o = 3 * (py * width + px);
      if (value !== value) {
        rgb[o] = rgb[o + 1] = rgb[o + 2] = 200;
        continue;
      }
      const k = Math.round(255 * map(value)) * 3;
      rgb[o] = table[k]; rgb[o + 1] = table[k + 1]; rgb[o + 2] = table[k + 2];
    }
  }
  return { png: encodePng(rgb, width, height), width, height };
}

/** The bin of a coordinate in ascending `edges`, or -1 outside them. */
function finder(edges) {
  const n = edges.length - 1, w = (edges[n] - edges[0]) / n;
  const uniform = edges.every((e, k) => Math.abs(e - edges[0] - k * w) <= 1e-6 * Math.abs(w));
  if (uniform) return (x) => (x >= edges[0] && x < edges[n] ? Math.min(n - 1, Math.floor((x - edges[0]) / w)) : -1);
  return (x) => {
    if (!(x >= edges[0] && x < edges[n])) return -1;
    let lo = 0, hi = n;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (x >= edges[mid]) lo = mid;
      else hi = mid;
    }
    return lo;
  };
}

// ---- PNG encoding (8-bit RGB, no interlace) ---------------------------------------------------

const CRC = Uint32Array.from({ length: 256 }, (_, n) => {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  return c >>> 0;
});

function crc32(bytes) {
  let c = 0xffffffff;
  for (const b of bytes) c = CRC[(c ^ b) & 255] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

function chunk(type, data) {
  const out = Buffer.alloc(12 + data.length);
  out.writeUInt32BE(data.length, 0);
  out.write(type, 4, 'latin1');
  Buffer.from(data).copy(out, 8);
  out.writeUInt32BE(crc32(out.subarray(4, 8 + data.length)), 8 + data.length);
  return out;
}

export function encodePng(rgb, width, height) {
  const header = Buffer.alloc(13);
  header.writeUInt32BE(width, 0);
  header.writeUInt32BE(height, 4);
  header.set([8, 2, 0, 0, 0], 8); // 8 bits, RGB
  const raw = Buffer.alloc((3 * width + 1) * height);
  for (let y = 0; y < height; y++) Buffer.from(rgb.buffer, rgb.byteOffset + 3 * width * y, 3 * width).copy(raw, y * (3 * width + 1) + 1);
  return new Uint8Array(Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk('IHDR', header), chunk('IDAT', deflateSync(raw)), chunk('IEND', Buffer.alloc(0)),
  ]));
}
