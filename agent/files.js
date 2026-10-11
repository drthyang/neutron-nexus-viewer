// Files the tools read and write: paths as the client gives them, .npy arrays, CSV text, and the rule
// that an existing file is replaced only when the call says so.

import { existsSync, mkdirSync, writeFileSync } from 'node:fs';
import { homedir } from 'node:os';
import { dirname, isAbsolute, resolve } from 'node:path';
import { gzipSync } from 'node:zlib';

import { ToolError } from './tool.js';

/** An absolute path: `~` is the home directory, a relative path starts from the server's working directory. */
export function resolvePath(path) {
  const p = String(path).trim();
  if (!p) throw new ToolError('Give a file path.');
  const home = p === '~' || p.startsWith('~/') ? homedir() + p.slice(1) : p;
  return isAbsolute(home) ? home : resolve(process.cwd(), home);
}

/**
 * Write `data` (bytes or text) to `path`, making its folder; `.gz` paths are gzipped. An existing file
 * is replaced only with `overwrite`. Returns { path, bytes }.
 */
export function writeOutput(path, data, overwrite = false) {
  const file = resolvePath(path);
  if (existsSync(file) && !overwrite) {
    throw new ToolError(`${file} exists. Pass overwrite: true to replace it, or another output_path.`);
  }
  let bytes = typeof data === 'string' ? Buffer.from(data, 'utf8') : Buffer.from(data.buffer, data.byteOffset, data.byteLength);
  if (/\.gz$/i.test(file)) bytes = gzipSync(bytes);
  mkdirSync(dirname(file), { recursive: true });
  writeFileSync(file, bytes);
  return { path: file, bytes: bytes.length };
}

const NPY_TYPES = new Map([
  [Float32Array, '<f4'], [Float64Array, '<f8'], [Uint8Array, '|u1'], [Uint16Array, '<u2'], [Uint32Array, '<u4'], [Int32Array, '<i4'],
]);

/** A NumPy .npy (version 1.0, C order) holding `array` with `shape`. */
export function npyBytes(array, shape) {
  const descr = NPY_TYPES.get(array.constructor);
  if (!descr) throw new Error(`No .npy type for ${array.constructor.name}.`);
  const dims = shape.length === 1 ? `${shape[0]},` : shape.join(', ');
  const dict = `{'descr': '${descr}', 'fortran_order': False, 'shape': (${dims}), }`;
  const padded = dict + ' '.repeat(63 - ((10 + dict.length) % 64)) + '\n';
  const out = new Uint8Array(10 + padded.length + array.byteLength);
  out.set([0x93, ...Buffer.from('NUMPY'), 1, 0, padded.length & 255, padded.length >> 8]);
  out.set(Buffer.from(padded, 'latin1'), 10);
  out.set(new Uint8Array(array.buffer, array.byteOffset, array.byteLength), 10 + padded.length);
  return out;
}

/** A number for a text file: up to 8 significant figures, NaN as `nan`. */
export const num = (x) => (Number.isFinite(x) ? String(Number(x.toPrecision(8))) : 'nan');
