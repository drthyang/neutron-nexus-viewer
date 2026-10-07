// Rigaku Oxford Diffraction (CrysAlisPro) experiment files: the binary frame
// format ("OD SAPPHIRE" header, TY6 compression, *.rod_img) and the text files
// that describe the instrument model and orientation (*.par, expinfo/*.ini).
//
// Header byte offsets follow FabIO (fabio/OXDimage.py) and dxtbx
// (format/FormatROD.py); the TY6 decoder implements the algorithm of dxtbx's
// decode_TY6_oneline. The decoder reproduces the min, max, mean and standard
// deviation that the instrument software stores in every frame header (checked
// on 3,638 HyPix-3000 frames; see docs/METHOD.md).

const ASCII = 256;
const text = (bytes) => String.fromCharCode(...bytes);

/** Parse the header of a frame (an ArrayBuffer or Uint8Array of at least the header size). */
export function parseRodHeader(buffer) {
  const u8 = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  const lines = text(u8.subarray(0, ASCII)).split('\r\n');
  const version = lines[0].trim();
  if (!/^OD /.test(version)) throw new Error('not a Rigaku Oxford Diffraction frame (no "OD" header)');
  const compression = (lines[1].split('=')[1] ?? '').slice(0, 3);
  const nx = parseInt(lines[2].slice(3, 7), 10), ny = parseInt(lines[2].slice(11, 15), 10);
  const l3 = lines[3];
  const nheader = parseInt(l3.slice(8, 15), 10);
  const sizes = {
    general: parseInt(l3.slice(19, 26), 10), special: parseInt(l3.slice(30, 37), 10), km4: parseInt(l3.slice(41, 48), 10),
    statistics: parseInt(l3.slice(52, 59), 10), history: parseInt(l3.slice(63), 10),
  };
  const time = (lines[5] ?? '').slice(5, 29).trim();
  const gen = ASCII, spe = gen + sizes.general, km4 = spe + sizes.special, sta = km4 + sizes.km4;
  const i16 = (o) => dv.getInt16(o, true), u16 = (o) => dv.getUint16(o, true);
  const i32 = (o) => dv.getInt32(o, true), u32 = (o) => dv.getUint32(o, true), f64 = (o) => dv.getFloat64(o, true);
  const npix = u32(gen + 36);
  if (npix !== nx * ny) throw new Error(`frame header: ${npix} pixels for ${nx} × ${ny}`);
  const steps = (o) => Array.from({ length: 10 }, (_, k) => i32(o + 4 * k));
  const stepDeg = Array.from({ length: 10 }, (_, k) => f64(km4 + 368 + 8 * k) * 180 / Math.PI);
  const deg = (s) => s.map((v, k) => v * stepDeg[k]);
  const start = deg(steps(km4 + 284)), end = deg(steps(km4 + 324));
  const moving = [0, 1, 2, 3].filter((k) => start[k] !== end[k]);
  return {
    version, compression, nx, ny, nheader, sizes, time,
    binning: [u16(gen), u16(gen + 2)],
    gain: f64(spe + 56),
    overflowFlag: i16(spe + 464),
    overflowThreshold: i32(spe + 472),
    exposure: f64(spe + 480),
    detectorType: i32(spe + 548),
    pixelMM: [f64(spe + 568), f64(spe + 576)],
    // Goniometer axes in header order: omega, theta, kappa, phi (+ six spare axes).
    start, end, stepDeg,
    zeroCorr: deg(steps(km4 + 512)),
    scanAxis: moving.length === 1 ? moving[0] : -1,
    beamRot: [f64(km4 + 552), f64(km4 + 560)],
    wavelengths: { alpha1: f64(km4 + 568), alpha2: f64(km4 + 576), alpha12: f64(km4 + 584), beta: f64(km4 + 592) },
    detRot: [f64(km4 + 640), f64(km4 + 648), f64(km4 + 656)],
    origin: [f64(km4 + 664), f64(km4 + 672)],
    alpha: f64(km4 + 680), beta: f64(km4 + 688),
    distance: f64(km4 + 712),
    stat: { min: i32(sta), max: i32(sta + 4), mean: f64(sta + 24), std: Math.sqrt(f64(sta + 32)) },
  };
}

/**
 * Decode the TY6-compressed pixels of a frame into `out` (Int32Array, ny*nx,
 * row-major: index = row * nx + column). `buffer` holds the whole file.
 */
export function decodeTY6(buffer, header, out = new Int32Array(header.nx * header.ny)) {
  const u8 = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  if (header.compression !== 'TY6') throw new Error(`compression ${header.compression} is not supported (TY6 only)`);
  const { nx, ny, nheader } = header;
  const nbytes = dv.getInt32(nheader, true);
  const base = nheader + 4, offsetsAt = base + nbytes;
  if (offsetsAt + 4 * ny > u8.length) throw new Error('truncated frame');
  const nblock = Math.floor((nx - 1) / 16), nrest = (nx - 1) % 16;
  const i16 = (p) => dv.getInt16(base + p, true), i32 = (p) => dv.getInt32(base + p, true);
  for (let row = 0; row < ny; row++) {
    let ip = dv.getUint32(offsetsAt + 4 * row, true);
    const o = row * nx;
    let op = o;
    const first = u8[base + ip++];
    if (first < 254) out[op] = first - 127;
    else if (first === 255) { out[op] = i32(ip); ip += 4; }
    else { out[op] = i16(ip); ip += 2; }
    op++;
    for (let k = 0; k < nblock; k++) {
      const bt = u8[base + ip++];
      for (let half = 0; half < 2; half++) {
        const nbit = (bt >> (4 * half)) & 15;
        if (nbit === 0) {
          for (let j = 0; j < 8; j++) out[op++] = 0;
          continue;
        }
        const zero = nbit > 1 ? (1 << (nbit - 1)) - 1 : 0;
        // 8 values of nbit bits packed little-endian in nbit bytes (nbit <= 8, so at most
        // 64 bits): the low 32 bits in lo, the rest in hi, read with unsigned shifts.
        let lo = 0, hi = 0;
        const p = base + ip;
        for (let j = 0; j < nbit && j < 4; j++) lo |= u8[p + j] << (8 * j);
        for (let j = 4; j < nbit; j++) hi |= u8[p + j] << (8 * (j - 4));
        ip += nbit;
        const mask = (1 << nbit) - 1;
        for (let j = 0, shift = 0; j < 8; j++, shift += nbit) {
          let field;
          if (shift + nbit <= 32) field = (lo >>> shift) & mask;
          else if (shift >= 32) field = (hi >>> (shift - 32)) & mask;
          else field = ((lo >>> shift) | (hi << (32 - shift))) & mask;
          out[op++] = field - zero;
        }
      }
      for (let i = op - 16; i < op; i++) {
        let d = out[i];
        if (d >= 127) {
          if (d >= 128) { d = i32(ip); ip += 4; } else { d = i16(ip); ip += 2; }
        }
        out[i] = d + out[i - 1];
      }
    }
    for (let i = 0; i < nrest; i++) {
      const px = u8[base + ip++];
      if (px < 254) out[op] = out[op - 1] + px - 127;
      else if (px === 255) { out[op] = out[op - 1] + i32(ip); ip += 4; }
      else { out[op] = out[op - 1] + i16(ip); ip += 2; }
      op++;
    }
    if (op !== o + nx) throw new Error(`TY6 row ${row}: decoded ${op - o} of ${nx} pixels`);
  }
  return out;
}

/** min, max, mean and population standard deviation of decoded pixels. */
export function frameStats(img) {
  let min = Infinity, max = -Infinity, s = 0, s2 = 0;
  for (let i = 0; i < img.length; i++) {
    const v = img[i];
    if (v < min) min = v;
    if (v > max) max = v;
    s += v;
    s2 += v * v;
  }
  const mean = s / img.length;
  return { min, max, mean, std: Math.sqrt(Math.max(0, s2 / img.length - mean * mean)) };
}

/** Does the decoded image reproduce the statistics stored in its header? */
export function checkStats(header, img) {
  const s = frameStats(img), h = header.stat;
  const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-12);
  return s.min === h.min && s.max === h.max && rel(s.mean, h.mean) < 1e-9 && rel(s.std, h.std) < 1e-6;
}

// ---- Experiment folder -----------------------------------------------------

const FRAME = /^(pre_)?(.+)_(\d+)_(\d+)\.rod_img$/;

/**
 * Sort the files of a CrysAlisPro experiment folder (objects with `name` and anything
 * else, passed through) into the main-series frames by run, the pre-experiment frames,
 * and the metadata files used here, found by name from the frames' experiment stem:
 * <stem>.par, <stem>_cracker.par, <stem>.p4p, <stem>_crystal.ini, <stem>_datacoll.ini.
 */
export function scanExperiment(files) {
  const runs = new Map(), pre = new Map(), stems = new Map();
  for (const f of files) {
    const m = FRAME.exec(f.name);
    if (!m) continue;
    const into = m[1] ? pre : runs;
    const run = Number(m[3]), frame = Number(m[4]);
    if (!into.has(run)) into.set(run, []);
    into.get(run).push({ frame, file: f, name: f.name });
    if (!m[1]) stems.set(m[2], (stems.get(m[2]) ?? 0) + 1);
  }
  if (stems.size > 1) throw new Error(`frames of several experiments in one folder: ${[...stems.keys()].join(', ')}`);
  const stem = stems.keys().next().value ?? null;
  const meta = {};
  const want = { par: `${stem}.par`, crackerPar: `${stem}_cracker.par`, p4p: `${stem}.p4p`, crystalIni: `${stem}_crystal.ini`, datacollIni: `${stem}_datacoll.ini` };
  for (const f of files) for (const [k, name] of Object.entries(want)) if (f.name === name) meta[k] ??= f;
  for (const m of [runs, pre]) for (const list of m.values()) list.sort((a, b) => a.frame - b.frame);
  const sorted = new Map([...runs.entries()].sort((a, b) => a[0] - b[0]));
  return { stem, runs: sorted, pre, meta };
}

/** Parse a Windows .ini file into { section: { key: value } } (keys as written). */
export function parseIni(txt) {
  const out = {};
  let section = '';
  for (const raw of txt.split(/\r?\n|\r/)) {
    const line = raw.trim();
    if (!line || line.startsWith(';') || line.startsWith('#')) continue;
    const s = /^\[(.*)\]$/.exec(line);
    if (s) { section = s[1]; out[section] ??= {}; continue; }
    const eq = line.indexOf('=');
    if (eq > 0) (out[section] ??= {})[line.slice(0, eq).trim()] = line.slice(eq + 1).trim();
  }
  return out;
}

const nums = (s) => (s ?? '').trim().split(/\s+/).map(Number).filter(Number.isFinite);

/**
 * crystal.ini: the orientation matrices CrysAlisPro stores (row-major, in
 * wavelength-scaled units), the cells, the Laue class and the instrument model
 * of the final refinement ("Mother of all lattices").
 */
export function parseCrystalIni(txt) {
  const ini = parseIni(txt);
  const ub = (sec) => {
    const v = nums(ini[sec]?.matrix);
    return v.length >= 9 ? v.slice(0, 9) : null;
  };
  const cell = (sec) => {
    const v = nums(ini[sec]?.['constants plus vol']);
    return v.length >= 6 && v[0] > 0 ? v.slice(0, 6) : null;
  };
  // InstrumentModel binary_0 is a hex dump of doubles; slots 27, 28 and 31 hold
  // the detector zero (x, y, pixels) and the distance (mm), 24-25 the detector rotations.
  const model = (hex) => {
    if (!hex) return null;
    const bytes = new Uint8Array(Math.floor(hex.length / 16) * 8);
    for (let i = 0; i < bytes.length; i++) bytes[i] = parseInt(hex.substr(2 * i, 2), 16);
    const d = new Float64Array(bytes.buffer);
    if (d.length < 32) return null;
    return { origin: [d[27], d[28]], distance: d[31], detRot: [d[24], d[25]] };
  };
  const laue = (ini.Symmetry?.['laue class'] ?? '').replace(/"/g, '').trim() || null;
  return {
    gralUB: ub('Gral UB'), motherUB: ub('Mother of all UBs'),
    cell: cell('Lattice'), constrainedCell: cell('Constrained lattice'), motherCell: cell('Mother of all lattices'),
    motherModel: model(ini['Mother of all lattices']?.['InstrumentModel binary_0']),
    laue, latticeType: (ini.Symmetry?.['lattice type'] ?? '').replace(/"/g, '').trim() || null,
  };
}

/**
 * A CrysAlisPro .par file: the UB (wavelength-scaled), the monochromator and
 * the temperature. Lines start with a non-ASCII marker for comments.
 */
export function parsePar(txt) {
  const out = { ub: null, monochromator: null, temperature: null, alpha: null, beta: null };
  for (const raw of txt.split(/\r?\n|\r/)) {
    const line = raw.replace(/^[^\x20-\x7e]+/, '').trim();
    let m;
    if ((m = /^CRYSTALLOGRAPHY UB\s+(.*)$/.exec(line))) {
      const v = nums(m[1]);
      if (v.length >= 9) out.ub = v.slice(0, 9);
    } else if ((m = /MONOCHROMATOR DVALUE \(ANG\)\s+([-\d.Ee+]+)\s+MONOCHROMATOR THETA \(DEG\)\s+([-\d.Ee+]+)/.exec(line))) {
      out.monochromator = { ...(out.monochromator ?? {}), dvalue: Number(m[1]), theta: Number(m[2]) };
    } else if ((m = /^ROTATION MONOCHROMATOR\s+([-\d.Ee+]+)\s+(\S+)/.exec(line))) {
      out.monochromator = { ...(out.monochromator ?? {}), plane: m[2] };
    } else if ((m = /^EXPERIMENT TEMPERATURE\s+([-\d.Ee+]+)/.exec(line))) {
      out.temperature = Number(m[1]);
    } else if ((m = /ALPHA \(DEG\)\s+([-\d.Ee+]+)\s+BETA \(DEG\)\s+([-\d.Ee+]+)/.exec(line))) {
      out.alpha = Number(m[1]);
      out.beta = Number(m[2]);
    }
  }
  return out;
}

/** expinfo/*_datacoll.ini: sample temperature range and timing, when recorded. */
export function parseDatacoll(txt) {
  const ini = parseIni(txt);
  const t = nums(ini['Sample T in K']?.['Sample T in K min max']);
  const unq = (s) => (s ?? '').replace(/"/g, '');
  return {
    temperature: t.length >= 2 ? t.slice(0, 2) : null,
    start: unq(ini.Date?.['Start time']) || null,
    end: unq(ini.Date?.['Successful end time']) || null,
  };
}
