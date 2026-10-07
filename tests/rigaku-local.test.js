// Checks against real CrysAlisPro data, which is not part of the repository: set
// RIGAKU_DIR to an experiment folder (with frames/*.rod_img) to run them, e.g.
//   RIGAKU_DIR=/path/to/experiment node --test tests/rigaku-local.test.js
// Every frame must decode to the min, max, mean and standard deviation stored in its header.

import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { checkStats, decodeTY6, parseRodHeader } from '../js/rigaku-format.js';

const dir = process.env.RIGAKU_DIR;

test('every frame of RIGAKU_DIR matches its header statistics', { skip: !dir && 'RIGAKU_DIR not set' }, () => {
  const frames = join(dir, 'frames');
  const names = readdirSync(frames).filter((n) => n.endsWith('.rod_img'));
  assert.ok(names.length > 0, 'no frames found');
  const bad = [];
  let out = null;
  for (const name of names) {
    const buf = readFileSync(join(frames, name));
    const h = parseRodHeader(buf);
    out = out?.length === h.nx * h.ny ? out : new Int32Array(h.nx * h.ny);
    decodeTY6(buf, h, out);
    if (!checkStats(h, out)) bad.push(name);
  }
  assert.deepEqual(bad, [], `${bad.length} of ${names.length} frames differ from their header statistics`);
});
