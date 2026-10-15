// Mantid MDHistoWorkspace files as Mantid's SaveMD (version 2) writes them, with h5wasm: the
// Rigaku reduction's and the reindexing's results, and the export for NEBULA3D.

import { cellFromUB, parseBasis } from './nexus.js';

/** The HDF5 inside h5wasm 0.10.3 (the version its error messages give; tests/rigaku.test.js keeps it in step). */
export const H5WASM_HDF5_VERSION = '2.0.0';

// Verbatim from a SaveMD file (CORELLI, Mantid 6), as NEBULA3D writes them: Sample.loadNexus
// parses the shape, and the parameter-map description is Mantid's own boilerplate.
const SHAPE_XML = '<type name="userShape">   <goniometer a11="1.000000" a12="0.000000" a13="0.000000" a21="0.000000" '
  + 'a22="1.000000" a23="0.000000" a31="0.000000" a32="0.000000" a33="1.000000"/></type>';
const PARAMETER_MAP_DESCRIPTION = 'A string representation of the parameter map. The format is either: '
  + '|detID:id-value;param-type;param-name;param-value| for a detector or  '
  + '|comp-name;param-type;param-name;param-value| for other components.';

const isAscii = (text) => /^[\x00-\x7f]*$/.test(text);

/** The local time as Mantid's file_time gives it: 2026-06-05T14:06:47-06:00. */
function fileTime(date) {
  const offset = -date.getTimezoneOffset(), two = (x) => String(Math.floor(Math.abs(x))).padStart(2, '0');
  const local = new Date(date.getTime() + offset * 60000).toISOString().slice(0, 19);
  return `${local}${offset < 0 ? '-' : '+'}${two(offset / 60)}:${two(offset % 60)}`;
}

const det3 = (m) => m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

/**
 * Write an MDHistoWorkspace as SaveMD v2 lays it out (checked against a file written by
 * Mantid 6.16 and reloaded with LoadMD): /MDHistoWorkspace/{data/{signal, errors_squared,
 * mask, num_events, D0, D1, D2}, experiment0/{instrument, logs, sample/oriented_lattice}}.
 *
 * `shape` is the storage shape (n_D2, n_D1, n_D0), so D0 is the fastest axis, and `labels`
 * name D0, D1, D2 in r.l.u. (by default H, K, L: storage (L, K, H)); `edges` are their bin
 * edges. `signal`, `errors2` (σ²), `events` (num_events) and `mask` (1 = masked; none masked
 * by default) are arrays in storage order, or functions (i, j) giving the slabs i to j along
 * the first axis. They are written one slab at a time and converted on the way, so no copy
 * of a whole volume in float64, or in the WASM heap, is ever made.
 *
 * `ub` is the UB matrix without 2π (Mantid's convention), row-major: the oriented lattice
 * when it is right-handed (Mantid refuses others) with `cell` (a, b, c, α, β, γ; from the UB
 * when omitted); without one, num_oriented_lattice is 0. `logs` become run logs, `fileName`
 * the root file_name attribute (the name the file is saved under), `gzip` the volumes'
 * compression level. `nebula3d` adds NEBULA3D's group of what Mantid does not know (`attrs`,
 * the exact bin `centers` of D0, D1, D2 and the `ub` with 2π), as NEBULA3D writes it, so
 * that NEBULA3D loads the file losslessly.
 *
 * The HDF5 types are Mantid's, since its loader is strict about them: the oldest file
 * format (superblock 0, as Mantid's files have), text attributes as variable-length ASCII
 * and text datasets as fixed-length null-terminated ASCII of shape (1,). Text outside ASCII
 * (a file name in a log, say) is labelled UTF-8 instead, so that every reader decodes it;
 * Mantid itself writes none. The volumes are gzip-compressed without SaveMD's byte shuffle,
 * which h5wasm cannot add, and `title` is the viewer's own: SaveMD writes none.
 */
export function writeMantidMD(h5wasm, path, {
  shape, edges, labels = ['[H,0,0]', '[0,K,0]', '[0,0,L]'], signal, errors2, events, mask = null,
  ub = null, cell = null, logs = {}, title = '', fileName = path.split('/').pop(), nebula3d = null, gzip = 4,
}) {
  const f = new h5wasm.File(path, 'w', { libver: ['earliest', 'latest'] });
  // h5wasm's dtype 'A' is ASCII (null-terminated), 'S' UTF-8; without a length, variable-length.
  const text = (obj, name, value) => obj.create_attribute(name, String(value), null, isAscii(String(value)) ? 'A' : 'S');
  const str = (g, name, value) => {
    const v = String(value) || ' ';
    g.create_dataset({ name, data: [v], shape: [1], dtype: `${isAscii(v) ? 'A' : 'S'}${new TextEncoder().encode(v).length}` });
  };
  const basis = labels.map((label) => parseBasis(label)?.vec);
  if (basis.some((v) => !v)) throw new Error(`Dimension labels ${labels.join(', ')} are not all r.l.u. directions such as [H,0,0].`);
  const lattice = ub && det3(ub) > 0 ? ub : null;
  const now = new Date();
  try {
    text(f, 'NX_class', 'NXroot');
    text(f, 'NeXus_version', '4.4.3');
    text(f, 'HDF5_Version', H5WASM_HDF5_VERSION);
    text(f, 'file_name', fileName);
    text(f, 'file_time', fileTime(now));
    const root = f.create_group('MDHistoWorkspace');
    text(root, 'NX_class', 'NXentry');
    text(root, 'QConvention', 'Crystallography');
    root.create_attribute('SaveMDVersion', 2, null, '<i4');
    root.create_dataset({ name: 'coordinate_system', data: Uint32Array.of(3) });
    root.create_dataset({ name: 'visual_normalization', data: Uint32Array.of(0) });
    if (title) str(root, 'title', title);
    const data = root.create_group('data');
    text(data, 'NX_class', 'NXdata');
    labels.forEach((long, i) => {
      const d = data.create_dataset({ name: `D${i}`, data: Float64Array.from(edges[i]) });
      text(d, 'frame', 'HKL');
      text(d, 'long_name', long);
      text(d, 'units', 'r.l.u.');
    });
    const plane = shape[1] * shape[2], chunks = [1, shape[1], shape[2]];
    const slabs = (x) => (typeof x === 'function' ? x : (i, j) => x.subarray(i * plane, j * plane));
    // One chunk (an axis-0 slab) per write: created from the first, resized, then filled. h5wasm's
    // dtypes are Python struct letters, whose size is the letter's: '<d' is float64 ('<f8' would
    // be float32), '|b' int8.
    const big = (name, dtype, source) => {
      const slab = slabs(source);
      const ds = data.create_dataset({ name, data: slab(0, 1), shape: chunks, maxshape: shape, chunks, dtype, compression: 'gzip', compression_opts: gzip });
      if (shape[0] > 1) ds.resize(shape);
      for (let i = 1; i < shape[0]; i++) ds.write_slice([[i, i + 1], [0, shape[1]], [0, shape[2]]], slab(i, i + 1));
      return ds;
    };
    const sig = big('signal', '<d', signal);
    text(sig, 'axes', 'D2:D1:D0');
    sig.create_attribute('signal', 1, null, '<i4');
    big('errors_squared', '<d', errors2);
    big('mask', '|b', mask ?? ((i, j) => new Int8Array((j - i) * plane)));
    big('num_events', '<d', events);
    const exp = root.create_group('experiment0');
    text(exp, 'NX_class', 'NXgroup');
    exp.create_attribute('version', 1, null, '<i4');
    const inst = exp.create_group('instrument');
    text(inst, 'NX_class', 'NXinstrument');
    inst.create_attribute('version', 1, null, '<i4');
    str(inst, 'name', ' ');
    const pmap = inst.create_group('instrument_parameter_map');
    text(pmap, 'NX_class', 'NXnote');
    pmap.create_attribute('version', 1, null, '<i4');
    // The date as Mantid writes it: UTC to the nanosecond, without a zone.
    const date = now.toISOString().slice(0, -1) + '000000';
    for (const [k, v] of [['author', ' '], ['data', ' '], ['date', date], ['description', PARAMETER_MAP_DESCRIPTION], ['type', 'text/plain']]) str(pmap, k, v);
    const xml = inst.create_group('instrument_xml');
    text(xml, 'NX_class', 'NXnote');
    for (const [k, v] of [['data', ' '], ['description', 'XML contents of the instrument IDF file.'], ['type', 'text/xml']]) str(xml, k, v);
    const lg = exp.create_group('logs');
    text(lg, 'NX_class', 'NXgroup');
    lg.create_attribute('version', 1, null, '<i4');
    const log = (name, value, units = ' ') => {
      const g = lg.create_group(name);
      text(g, 'NX_class', 'NXlog');
      const d = typeof value === 'string' ? (str(g, 'value', value), g.get('value'))
        : g.create_dataset({ name: 'value', data: Float64Array.from([].concat(value)) });
      text(d, 'units', units);
    };
    // Column j is the (h, k, l) direction of Dj.
    log('W_MATRIX', [0, 1, 2].flatMap((i) => basis.map((v) => v[i])));
    const gon = lg.create_group('goniometer');
    text(gon, 'NX_class', 'NXpositioner');
    gon.create_attribute('version', 1, null, '<i4');
    gon.create_dataset({ name: 'num_axes', data: Int32Array.of(0) });
    gon.create_dataset({ name: 'rotation_matrix', data: Float64Array.of(1, 0, 0, 0, 1, 0, 0, 0, 1) });
    log('mdhisto_was_modified', '1');
    for (const [name, v] of Object.entries(logs)) {
      if (v && typeof v === 'object' && 'value' in v) log(name, v.value, v.units ?? ' ');
      else log(name, typeof v === 'number' ? v : String(v));
    }
    const sample = exp.create_group('sample');
    text(sample, 'NX_class', 'NXsample');
    text(sample, 'name', ' ');
    sample.create_attribute('name_empty', 1, null, '<i4');
    text(sample, 'shape_xml', SHAPE_XML);
    sample.create_attribute('version', 1, null, '<i4');
    for (const k of ['geom_height', 'geom_thickness', 'geom_width']) sample.create_dataset({ name: k, data: Float64Array.of(0) });
    sample.create_dataset({ name: 'geom_id', data: Int32Array.of(0) });
    const mat = sample.create_group('material');
    text(mat, 'NX_class', 'NXdata');
    text(mat, 'formulaStyle', 'empty');
    text(mat, 'name', ' ');
    mat.create_attribute('version', 2, null, '<i4');
    for (const [k, v] of [['number_density', 0], ['packing_fraction', 1], ['pressure', 0], ['temperature', 0]]) mat.create_dataset({ name: k, data: Float64Array.of(v) });
    sample.create_dataset({ name: 'num_oriented_lattice', data: Int32Array.of(lattice ? 1 : 0) });
    sample.create_dataset({ name: 'num_other_samples', data: Int32Array.of(0) });
    if (lattice) {
      const c = cell ?? (({ a, b, c: cc, alpha, beta, gamma }) => [a, b, cc, alpha, beta, gamma])(cellFromUB(lattice));
      const ol = sample.create_group('oriented_lattice');
      text(ol, 'NX_class', 'NXcrystal');
      ol.create_dataset({ name: 'cross_term', data: Int32Array.of(0) });
      ol.create_dataset({ name: 'maximum_order', data: Int32Array.of(0) });
      for (const k of ['modulated_hkl_error', 'modulated_orientation_matrix']) ol.create_dataset({ name: k, data: new Float64Array(9), shape: [3, 3] });
      ol.create_dataset({ name: 'orientation_matrix', data: Float64Array.from(lattice), shape: [3, 3] });
      ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].forEach((k, i) => {
        ol.create_dataset({ name: `unit_cell_${k}`, data: Float64Array.of(c[i]) });
        ol.create_dataset({ name: `unit_cell_${k}_error`, data: Float64Array.of(0) });
      });
    }
    if (nebula3d) {
      // NEBULA3D's own group (its format version 1): LoadMD opens only the groups it knows.
      const extra = root.create_group('nebula3d');
      text(extra, 'NX_class', 'NXcollection');
      extra.create_attribute('format_version', 1, null, '<i4');
      for (const [k, v] of Object.entries(nebula3d.attrs ?? {})) text(extra, k, v);
      nebula3d.centers.forEach((c, i) => extra.create_dataset({ name: `D${i}_centers`, data: Float64Array.from(c) }));
      if (nebula3d.ub) extra.create_dataset({ name: 'ub_matrix', data: Float64Array.from(nebula3d.ub), shape: [3, 3] });
    }
  } finally {
    f.close();
  }
}
