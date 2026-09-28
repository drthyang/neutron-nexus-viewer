import assert from 'node:assert/strict';
import test from 'node:test';

import { closeGroup, formatOp, indexMaps, metricChange, parseOp, parseOps, PRESETS } from '../js/symmetry.js';

const group = (name) => closeGroup(parseOps(PRESETS.find(([n]) => n === name)[1]));
const apply = (R, v) => [0, 1, 2].map((i) => R[3 * i] * v[0] + R[3 * i + 1] * v[1] + R[3 * i + 2] * v[2]);
const equivalents = (name, v) => new Set(group(name).map((R) => apply(R, v).join(',')));

test('Laue class orders', () => {
  const orders = {
    1: 1, '-1': 2, '2/m (b unique)': 4, '2/m (c unique)': 4, mmm: 8, '4/m': 8, '4/mmm': 16,
    '-3': 6, '-3m1': 12, '-31m': 12, '6/m': 12, '6/mmm': 24, 'm-3': 24, 'm-3m': 48,
  };
  for (const [name] of PRESETS) assert.equal(group(name).length, orders[name], name);
});

test('equivalent reflections', () => {
  assert.deepEqual(equivalents('6/mmm', [1, 0, 0]), new Set(['1,0,0', '0,1,0', '-1,1,0', '-1,0,0', '0,-1,0', '1,-1,0']));
  assert.equal(equivalents('6/mmm', [1, 2, 3]).size, 24);
  // -3m1 relates (h,k,l) to (k,h,-l); -31m relates it to (k,h,l).
  assert.ok(equivalents('-3m1', [1, 2, 3]).has('2,1,-3'));
  assert.ok(!equivalents('-3m1', [1, 2, 3]).has('2,1,3'));
  assert.ok(equivalents('-31m', [1, 2, 3]).has('2,1,3'));
  assert.equal(equivalents('m-3m', [1, 2, 3]).size, 48);
  assert.equal(equivalents('m-3', [1, 2, 3]).size, 24);
});

test('parsing and formatting', () => {
  assert.deepEqual(parseOp('h+k, -h, l'), [1, 1, 0, -1, 0, 0, 0, 0, 1]);
  assert.deepEqual(parseOp('-h+k,k,-l').slice(0, 3), [-1, 1, 0]);
  // Real-space x-y,x,z (6+) acts on reflections as its transpose.
  assert.deepEqual(parseOp('x-y,x,z'), parseOp('h+k,-h,l'));
  assert.equal(formatOp(parseOp('h+k,-h,l')), 'h+k,-h,l');
  assert.equal(formatOp(parseOp('-h+k,k,-l')), '-h+k,k,-l');
  assert.deepEqual(parseOps('None'), []);
  assert.deepEqual(parseOps('-h,-k,-l; k,h,l\nh,k,-l').length, 3);
});

test('invalid operations are rejected', () => {
  assert.throws(() => parseOp('h,k'), /three/);
  assert.throws(() => parseOp('h+1/2,k,l'), /translations/);
  assert.throws(() => parseOp('2h,k,l'), /determinant/);
  assert.throws(() => parseOp('h,y,l'), /mixes/);
  assert.throws(() => closeGroup([parseOp('h+k,k,l')]), /not a crystallographic point group/);
});

test('metric change flags operations that are not cell symmetries', () => {
  const hex = { a: 8, b: 8, c: 10, alpha: 90, beta: 90, gamma: 120 };
  assert.ok(metricChange(group('6/mmm'), hex) < 1e-12);
  assert.ok(metricChange(group('mmm'), hex) > 0.1);
  const cubic = { a: 5, b: 5, c: 5, alpha: 90, beta: 90, gamma: 90 };
  assert.ok(metricChange(group('m-3m'), cubic) < 1e-12);
});

test('index maps require a compatible grid', () => {
  const axis = (n, w, name) => ({
    label: name, basis: { vec: [...'HKL'].map((c) => (c === name ? 1 : 0)), letter: name },
    edges: Array.from({ length: n + 1 }, (_, i) => (i - n / 2) * w),
  });
  const good = [axis(5, 0.1, 'H'), axis(5, 0.1, 'K'), axis(3, 0.2, 'L')];
  const [six] = indexMaps([parseOp('h+k,-h,l')], good);
  // Index (i, j, m) with the origin at (2, 2, 1): h+k -> i+j-2, -h -> 4-i, l -> m.
  assert.deepEqual([...six.M], [1, 1, 0, -1, 0, 0, 0, 0, 1]);
  assert.deepEqual([...six.t], [-2, 4, 0]);
  assert.throws(() => indexMaps([parseOp('h+k,-h,l')], [axis(5, 0.1, 'H'), axis(5, 0.2, 'K'), axis(3, 0.2, 'L')]), /bin widths/);
  assert.throws(() => indexMaps([parseOp('h+k,-h,l')], [axis(4, 0.1, 'H'), axis(4, 0.1, 'K'), axis(3, 0.2, 'L')]), /centered at 0/);
  // Operations follow the axis basis: with axes [H,H,0], [-K,K,0], L the 2-fold
  // (k,h,l) is a mirror of the second axis.
  const rotated = [
    { label: 'HH0', basis: { vec: [1, 1, 0] }, edges: good[0].edges },
    { label: '-KK0', basis: { vec: [-1, 1, 0] }, edges: good[1].edges },
    good[2],
  ];
  const [swap] = indexMaps([parseOp('k,h,l')], rotated);
  assert.deepEqual([...swap.M], [1, 0, 0, 0, -1, 0, 0, 0, 1]);
});
