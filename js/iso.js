// Coarse volume and isosurface for the 3-D view.
//
// The volume is binned by an odd factor with blocks aligned to the origin, so
// symmetry operations that map the fine grid onto itself also map the coarse
// grid onto itself. Blocks pool their valid voxels (sum and count), orbits pool
// their unique blocks, and the isosurface of the pooled mean is extracted with
// naive surface nets (one vertex per boundary cell, one quad per crossed edge).

/** Coarse binning of display dims so no axis has more than maxBins blocks. */
export function coarseGrid(dims, maxBins) {
  const n = dims.map((d) => d.edges.length - 1);
  let factor = Math.max(1, Math.ceil(Math.max(...n) / maxBins));
  if (factor % 2 === 0) factor++;
  const axes = dims.map((d, i) => {
    const w = (d.edges[n[i]] - d.edges[0]) / n[i];
    const centered = Math.abs(d.edges[0] + d.edges[n[i]]) < 1e-4 * w;
    let start = 0, count = Math.ceil(n[i] / factor);
    if (centered && n[i] % 2) {
      // A block centered on the middle bin, which sits at the origin.
      const c = (n[i] - 1) / 2, half = (factor - 1) / 2, k = Math.ceil((c - half) / factor);
      start = c - half - k * factor;
      count = 2 * k + 1;
    } else if (centered) {
      // A block boundary on the origin, which is a bin edge.
      const k = Math.ceil(n[i] / 2 / factor);
      start = n[i] / 2 - k * factor;
      count = 2 * k;
    }
    const edges = Array.from({ length: count + 1 }, (_, m) => d.edges[0] + (start + m * factor) * w);
    const block = Int32Array.from({ length: n[i] }, (_, j) => Math.floor((j - start) / factor));
    return { start, count, block, dim: { ...d, edges } };
  });
  return { factor, shape: axes.map((a) => a.count), block: axes.map((a) => a.block), dims: axes.map((a) => a.dim) };
}

/** Sum and count valid (non-NaN) voxels per coarse block. Coarse arrays use display order (dim 0 fastest). */
export function binVolume(volume, shape, grid) {
  const [n0, n1, n2] = [shape[2], shape[1], shape[0]];
  const [c0, c1] = grid.shape, total = grid.shape[0] * grid.shape[1] * grid.shape[2];
  const [b0, b1, b2] = grid.block;
  const sums = new Float64Array(total), counts = new Uint32Array(total);
  for (let i2 = 0, f = 0; i2 < n2; i2++) {
    for (let i1 = 0; i1 < n1; i1++) {
      const base = (b2[i2] * c1 + b1[i1]) * c0;
      for (let i0 = 0; i0 < n0; i0++, f++) {
        const v = volume[f];
        if (v === v) { const k = base + b0[i0]; sums[k] += v; counts[k]++; }
      }
    }
  }
  return { sums, counts };
}

/**
 * Pooled mean per coarse voxel over its orbit's unique voxels (NaN where the
 * orbit has no data). `maps` are index maps on the coarse grid, identity first.
 */
export function orbitMean({ sums, counts }, shape, maps) {
  const [n0, n1, n2] = shape, mean = new Float32Array(sums.length);
  for (let i2 = 0, p = 0; i2 < n2; i2++) {
    for (let i1 = 0; i1 < n1; i1++) {
      for (let i0 = 0; i0 < n0; i0++, p++) {
        let s = 0, c = 0, stab = 0;
        for (const { M, t } of maps) {
          const q0 = M[0] * i0 + M[1] * i1 + M[2] * i2 + t[0];
          const q1 = M[3] * i0 + M[4] * i1 + M[5] * i2 + t[1];
          const q2 = M[6] * i0 + M[7] * i1 + M[8] * i2 + t[2];
          if (q0 === i0 && q1 === i1 && q2 === i2) stab++;
          if (q0 < 0 || q1 < 0 || q2 < 0 || q0 >= n0 || q1 >= n1 || q2 >= n2) continue;
          const q = (q2 * n1 + q1) * n0 + q0;
          s += sums[q];
          c += counts[q];
        }
        // Each unique image appears |stabilizer| times; the ratio is unaffected.
        mean[p] = c ? s / c : NaN;
      }
    }
  }
  return mean;
}

const CORNERS = Array.from({ length: 8 }, (_, c) => [c & 1, (c >> 1) & 1, (c >> 2) & 1]);
const EDGES = [];
for (let c = 0; c < 8; c++) for (const bit of [1, 2, 4]) if (!(c & bit)) EDGES.push([c, c | bit]);

/**
 * Surface nets on a scalar field in display order (dim 0 fastest). NaN counts
 * as below the level, so surfaces close where there is no data. Returns vertex
 * positions in grid-index coordinates and triangle indices.
 */
export function surfaceNets(field, shape, level) {
  const [n0, n1, n2] = shape;
  const value = (i) => { const v = field[i]; return v === v ? v : -Infinity; };
  const cells = new Int32Array(Math.max(0, (n0 - 1) * (n1 - 1) * (n2 - 1))).fill(-1);
  const cellIndex = (i, j, k) => (k * (n1 - 1) + j) * (n0 - 1) + i;
  const positions = [], indices = [];
  const corner = new Float64Array(8);

  for (let k = 0; k < n2 - 1; k++) {
    for (let j = 0; j < n1 - 1; j++) {
      for (let i = 0; i < n0 - 1; i++) {
        let mask = 0;
        for (let c = 0; c < 8; c++) {
          const [dx, dy, dz] = CORNERS[c];
          corner[c] = value(((k + dz) * n1 + j + dy) * n0 + i + dx);
          if (corner[c] > level) mask |= 1 << c;
        }
        if (mask === 0 || mask === 255) continue;
        let x = 0, y = 0, z = 0, crossings = 0;
        for (const [a, b] of EDGES) {
          if (((mask >> a) & 1) === ((mask >> b) & 1)) continue;
          const va = corner[a], vb = corner[b];
          const t = va === -Infinity ? 1 : vb === -Infinity ? 0 : (level - va) / (vb - va);
          x += CORNERS[a][0] + t * (CORNERS[b][0] - CORNERS[a][0]);
          y += CORNERS[a][1] + t * (CORNERS[b][1] - CORNERS[a][1]);
          z += CORNERS[a][2] + t * (CORNERS[b][2] - CORNERS[a][2]);
          crossings++;
        }
        cells[cellIndex(i, j, k)] = positions.length / 3;
        positions.push(i + x / crossings, j + y / crossings, k + z / crossings);
      }
    }
  }

  // One quad per grid edge whose ends straddle the level, joining the four
  // cells around that edge; winding follows which end is inside.
  const n = [n0, n1, n2];
  for (let a = 0; a < 3; a++) {
    const b = (a + 1) % 3, c = (a + 2) % 3;
    const step = [1, n0, n0 * n1];
    const p = [0, 0, 0];
    for (p[2] = 0; p[2] < n2; p[2]++) {
      for (p[1] = 0; p[1] < n1; p[1]++) {
        for (p[0] = 0; p[0] < n0; p[0]++) {
          if (p[a] >= n[a] - 1 || p[b] < 1 || p[c] < 1 || p[b] >= n[b] - 1 || p[c] >= n[c] - 1) continue;
          const node = p[0] + step[1] * p[1] + step[2] * p[2];
          const inside = value(node) > level;
          if (inside === value(node + step[a]) > level) continue;
          const cell = (db, dc) => {
            const q = [...p];
            q[b] -= db;
            q[c] -= dc;
            return cells[cellIndex(q[0], q[1], q[2])];
          };
          const v00 = cell(1, 1), v10 = cell(0, 1), v11 = cell(0, 0), v01 = cell(1, 0);
          if (inside) indices.push(v00, v10, v11, v00, v11, v01);
          else indices.push(v00, v11, v10, v00, v01, v11);
        }
      }
    }
  }
  return { positions: Float32Array.from(positions), indices: Uint32Array.from(indices) };
}
