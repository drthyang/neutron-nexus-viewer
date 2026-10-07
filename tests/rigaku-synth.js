// Test helpers: write synthetic Rigaku Oxford Diffraction frames ("OD SAPPHIRE 4.0",
// TY6) and simulate a rotation experiment with a known geometry, for rigaku.test.js.
// The TY6 encoder produces the byte stream that js/rigaku-format.js decodes (the
// decoder itself is checked against real HyPix-3000 frames; see docs/METHOD.md).

import { M_CRYSALIS, bMatrix, mul, mulv, predict, prepare, rot, transpose } from '../js/rigaku-geometry.js';

/** TY6-encode a ny x nx Int32 image: { data (Uint8Array), offsets (Uint32Array) }. */
export function encodeTY6(img, nx, ny) {
  const out = [], offsets = new Uint32Array(ny);
  const i16 = (v) => out.push(v & 255, (v >> 8) & 255);
  const i32 = (v) => out.push(v & 255, (v >> 8) & 255, (v >> 16) & 255, (v >>> 24) & 255);
  const single = (d) => {
    if (d >= -127 && d <= 126) out.push(d + 127);
    else if (d >= -32768 && d <= 32767) { out.push(254); i16(d); } else { out.push(255); i32(d); }
  };
  for (let y = 0; y < ny; y++) {
    offsets[y] = out.length;
    const row = img.subarray(y * nx, (y + 1) * nx);
    single(row[0]);
    const d = [];
    for (let x = 1; x < nx; x++) d.push(row[x] - row[x - 1]);
    const nblock = Math.floor((nx - 1) / 16);
    for (let b = 0; b < nblock; b++) {
      const block = d.slice(16 * b, 16 * b + 16);
      const halves = [block.slice(0, 8), block.slice(8)], nbits = [], fields = [], extra = [];
      for (const h of halves) {
        let nbit = 8;
        for (let n = 0; n < 8; n++) {
          const zero = n > 1 ? (1 << (n - 1)) - 1 : 0, hi = (1 << n) - 1 - zero;
          if (h.every((v) => v >= -zero && v <= hi)) { nbit = n; break; }
        }
        const zero = nbit > 1 ? (1 << (nbit - 1)) - 1 : 0;
        nbits.push(nbit);
        fields.push(h.map((v) => {
          if (nbit < 8 || (v >= -127 && v <= 126)) return v + zero;
          extra.push(v);
          return v >= -32768 && v <= 32767 ? 254 : 255;
        }));
      }
      out.push(nbits[0] | (nbits[1] << 4));
      halves.forEach((_, k) => {
        const nbit = nbits[k];
        let bits = 0n;
        fields[k].forEach((f, j) => { bits |= BigInt(f) << BigInt(nbit * j); });
        for (let j = 0; j < nbit; j++) out.push(Number((bits >> BigInt(8 * j)) & 255n));
      });
      for (const v of extra) {
        if (v >= -32768 && v <= 32767) i16(v); else i32(v);
      }
    }
    for (const v of d.slice(16 * nblock)) single(v);
  }
  return { data: Uint8Array.from(out), offsets };
}

const pad = (v, n) => String(v).padStart(n, ' ');

/**
 * A complete frame file. `h`: { nx, ny, start: [omega, theta, kappa, phi], end, zeroCorr,
 * exposure, pixelMM, origin: [x, y], distance, detRot: [d1, d2, d3], alpha, beta,
 * wavelengths: [a1, a2, a12, beta], overflowThreshold }.
 */
export function writeFrame(img, h) {
  const { nx, ny } = h;
  const NHEADER = 6576, sizes = [512, 768, 1024, 512, 2048];
  const head = new Uint8Array(NHEADER), dv = new DataView(head.buffer);
  const ascii = ['OD SAPPHIRE  4.0', 'COMPRESSION=TY6(  7.6)', `NX=${pad(nx, 4)} NY=${pad(ny, 4)} OI=${pad(0, 7)} OL=${pad(0, 7)} `,
    `NHEADER=${pad(NHEADER, 7)} NG=${pad(sizes[0], 7)} NS=${pad(sizes[1], 7)} NK=${pad(sizes[2], 7)} NS=${pad(sizes[3], 7)} NH=${pad(sizes[4], 7)}`,
    `NSUPPLEMENT=${pad(0, 7)}`, 'TIME=Thu Jan 01 00:00:00 2026'].join('\r\n') + '\r\n';
  head.set(new TextEncoder().encode(ascii.padEnd(255, ' ')).subarray(0, 255));
  head[255] = 0x1a;
  const gen = 256, spe = gen + 512, km4 = spe + 768, sta = km4 + 1024;
  dv.setUint16(gen, 1, true); dv.setUint16(gen + 2, 1, true);
  dv.setUint32(gen + 36, nx * ny, true);
  dv.setFloat64(spe + 56, 1, true);
  dv.setInt32(spe + 472, h.overflowThreshold ?? 1e6, true);
  dv.setFloat64(spe + 480, h.exposure, true);
  dv.setInt32(spe + 548, 7, true);
  dv.setFloat64(spe + 568, h.pixelMM, true); dv.setFloat64(spe + 576, h.pixelMM, true);
  const step = 1e-4; // degrees per motor step
  for (let k = 0; k < 4; k++) {
    dv.setInt32(km4 + 284 + 4 * k, Math.round(h.start[k] / step), true);
    dv.setInt32(km4 + 324 + 4 * k, Math.round(h.end[k] / step), true);
    dv.setInt32(km4 + 512 + 4 * k, Math.round((h.zeroCorr?.[k] ?? 0) / step), true);
  }
  for (let k = 0; k < 10; k++) dv.setFloat64(km4 + 368 + 8 * k, step * Math.PI / 180, true);
  h.wavelengths.forEach((w, k) => dv.setFloat64(km4 + 568 + 8 * k, w, true));
  h.detRot.forEach((v, k) => dv.setFloat64(km4 + 640 + 8 * k, v, true));
  dv.setFloat64(km4 + 664, h.origin[0], true); dv.setFloat64(km4 + 672, h.origin[1], true);
  dv.setFloat64(km4 + 680, h.alpha, true); dv.setFloat64(km4 + 688, h.beta ?? 0, true);
  dv.setFloat64(km4 + 712, h.distance, true);
  let min = Infinity, max = -Infinity, s = 0, s2 = 0;
  for (const v of img) { min = Math.min(min, v); max = Math.max(max, v); s += v; s2 += v * v; }
  const mean = s / img.length;
  dv.setInt32(sta, min, true); dv.setInt32(sta + 4, max, true);
  dv.setFloat64(sta + 24, mean, true); dv.setFloat64(sta + 32, s2 / img.length - mean * mean, true);
  const { data, offsets } = encodeTY6(img, nx, ny);
  const file = new Uint8Array(NHEADER + 4 + data.length + 4 * ny);
  file.set(head);
  new DataView(file.buffer).setInt32(NHEADER, data.length, true);
  file.set(data, NHEADER + 4);
  file.set(new Uint8Array(offsets.buffer), NHEADER + 4 + data.length);
  return file;
}

/** Deterministic pseudo-random numbers (mulberry32). */
export function rng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const erf = (x) => {
  // Abramowitz-Stegun 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
};

/**
 * Simulate an omega-scan experiment with a hexagonal crystal: Gaussian Bragg spots on a
 * flat background, frames as `{ run, frame, bytes }`. Returns the true model too.
 */
export function simulateExperiment({ nx = 160, ny = 100, runs = [[54, 0], [54, 90]], width = 1, start = -70, nFrames = 170, sigmaOmega = 0.6,
  cell = [4.0, 4.0, 5.0, 90, 90, 120], U = mul(rot([1, 0, 0], 23), mul(rot([0, 1, 0], 41), rot([0, 0, 1], 17))),
  geometry = {}, amplitude = 4000, background = 2, seed = 1 } = {}) {
  const g = {
    wavelength: 0.71073, ox: nx / 2 + 3.3, oy: ny / 2 - 2.1, distance: 45, pixelMM: 0.4, d1: 0.3, d2: -0.1,
    theta: 20, thetaOffset: 0, alpha: 90, omegaOffset: 90, kappaOffset: 0, phiOffset: 0, tx: 0, ty: 0, tz: 0, ...geometry,
  };
  const G = prepare(g);
  const ub = mul(U, bMatrix(...cell)).map((v) => v * g.wavelength);
  const MT = transpose(M_CRYSALIS);
  const spots = [];
  for (const [r, [kappa, phi]] of runs.entries()) {
    for (let h = -6; h <= 6; h++) for (let k = -6; k <= 6; k++) for (let l = -7; l <= 7; l++) {
      if (!h && !k && !l) continue;
      const xC = mulv(mul(MT, ub), [h, k, l]);
      const found = new Set();
      for (let near = start; near <= start + width * nFrames; near += 20) {
        const p = predict(G, xC, [0, g.theta, kappa, phi], 0, near);
        if (!p || p.angle < start || p.angle > start + width * nFrames) continue;
        const key = p.angle.toFixed(3);
        if (found.has(key)) continue;
        found.add(key);
        if (p.i > 2 && p.j > 2 && p.i < nx - 3 && p.j < ny - 3) spots.push({ run: r + 1, hkl: [h, k, l], ...p });
      }
    }
  }
  const rand = rng(seed), frames = [];
  for (const [r, [kappa, phi]] of runs.entries()) {
    for (let f = 0; f < nFrames; f++) {
      const a0 = start + f * width, a1 = a0 + width;
      const img = new Int32Array(nx * ny);
      for (let q = 0; q < img.length; q++) img[q] = Math.round(background + (rand() - 0.5) * 2);
      for (const s of spots) {
        if (s.run !== r + 1) continue;
        const wf = 0.5 * (erf((a1 - s.angle) / (sigmaOmega * Math.SQRT2)) - erf((a0 - s.angle) / (sigmaOmega * Math.SQRT2)));
        if (wf < 1e-4) continue;
        for (let y = Math.floor(s.j) - 3; y <= Math.ceil(s.j) + 3; y++) {
          for (let x = Math.floor(s.i) - 3; x <= Math.ceil(s.i) + 3; x++) {
            if (x < 0 || y < 0 || x >= nx || y >= ny) continue;
            img[y * nx + x] += Math.round(amplitude * wf * Math.exp(-((x - s.i) ** 2 + (y - s.j) ** 2) / (2 * 0.8 ** 2)) / (2 * Math.PI * 0.64));
          }
        }
      }
      const bytes = writeFrame(img, {
        nx, ny, start: [a0, g.theta, kappa, phi], end: [a1, g.theta, kappa, phi], zeroCorr: [g.omegaOffset, 0, 0, 0], exposure: 1,
        pixelMM: g.pixelMM, origin: [g.ox, g.oy], distance: g.distance, detRot: [g.d1, g.d2, 0], alpha: g.alpha,
        wavelengths: [0.7093, 0.71359, g.wavelength, 0.63229],
      });
      frames.push({ run: r + 1, frame: f + 1, bytes, img });
    }
  }
  return { frames, spots, model: { g, ub, cell } };
}
