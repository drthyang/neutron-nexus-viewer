// Geometry of a Rigaku Oxford Diffraction (CrysAlisPro) kappa goniometer with an
// area detector on a 2theta arm, and refinement of it against Bragg peaks.
//
// Laboratory frame (as dxtbx): +Z from the sample towards the source (the beam
// travels along -Z, s0 = (0, 0, -1)), +Y = -omega axis, X = Y x Z. Rotations are
// right-handed.
//   goniometer  x_lab = R(e_omega, w) R(e_kappa, k) R(e_phi, p) x_C,
//               e_omega = e_phi = (0, -1, 0), e_kappa = (0, -cos(alpha), sin(alpha)),
//               w = omega_header + omegaOffset, k = kappa_header + kappaOffset;
//               omegaOffset starts at the header's software zero correction (+90 deg
//               on the XtaLAB mini II).
//   detector    R_det = R((0,-1,0), theta + thetaOffset) R((-1,0,0), d2) R((0,0,1), d1);
//               pixel (i, j) centre -> R_det [(i - ox) px, (j - oy) px, -D].
//   reciprocal  x = s - s0 (|x| = 2 sin(theta_B) = lambda/d: "lambda units", no 2*pi).
//               The CrysAlis UB is in the same units, in its own frame e1 = +Z,
//               e2 = +X, e3 = +Y:  hkl = UB^-1 M x_C, M = [[0,0,1],[1,0,0],[0,1,0]].
// These conventions were established empirically on two XtaLAB mini II datasets:
// out of 4,608 and 6,144 combinations of axis permutations, offsets and rotation
// senses, only these index (nearly) all peaks of all runs (docs/METHOD.md).

export const M_CRYSALIS = [0, 0, 1, 1, 0, 0, 0, 1, 0];
export const DEG = Math.PI / 180;

// ---- 3x3 matrices, row-major arrays of 9 -----------------------------------

export function mul(a, b) {
  const c = new Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) c[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
  }
  return c;
}
export const mulv = (a, v) => [a[0] * v[0] + a[1] * v[1] + a[2] * v[2], a[3] * v[0] + a[4] * v[1] + a[5] * v[2],
  a[6] * v[0] + a[7] * v[1] + a[8] * v[2]];
export const transpose = (a) => [a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8]];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);

export function inv3(m) {
  const d = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
  if (!(Math.abs(d) > 1e-300)) throw new Error('singular matrix');
  return [(m[4] * m[8] - m[5] * m[7]) / d, (m[2] * m[7] - m[1] * m[8]) / d, (m[1] * m[5] - m[2] * m[4]) / d,
    (m[5] * m[6] - m[3] * m[8]) / d, (m[0] * m[8] - m[2] * m[6]) / d, (m[2] * m[3] - m[0] * m[5]) / d,
    (m[3] * m[7] - m[4] * m[6]) / d, (m[1] * m[6] - m[0] * m[7]) / d, (m[0] * m[4] - m[1] * m[3]) / d];
}

/** Right-handed rotation by `deg` degrees about `axis`. */
export function rot(axis, deg) {
  const n = norm(axis), [x, y, z] = axis.map((v) => v / n);
  const t = deg * DEG, c = Math.cos(t), s = Math.sin(t), C = 1 - c;
  return [c + x * x * C, x * y * C - z * s, x * z * C + y * s,
    y * x * C + z * s, c + y * y * C, y * z * C - x * s,
    z * x * C - y * s, z * y * C + x * s, c + z * z * C];
}

/** Small rotation (degrees) about X, then Y, then Z composed as Rx Ry Rz. */
export const smallRot = ([a, b, c]) => mul(mul(rot([1, 0, 0], a), rot([0, 1, 0], b)), rot([0, 0, 1], c));

// ---- Cells --------------------------------------------------------------------

/** Busing-Levy B matrix (1/Angstrom, no 2*pi) of a direct cell. */
export function bMatrix(a, b, c, al, be, ga) {
  const [ca, cb, cg] = [al, be, ga].map((x) => Math.cos(x * DEG));
  const [sa, sb, sg] = [al, be, ga].map((x) => Math.sin(x * DEG));
  const V = a * b * c * Math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg);
  const as = b * c * sa / V, bs = a * c * sb / V, cs = a * b * sg / V;
  const cal = (cb * cg - ca) / (sb * sg), cbe = (ca * cg - cb) / (sa * sg), cga = (ca * cb - cg) / (sa * sb);
  const sbe = Math.sqrt(1 - cbe * cbe), sga = Math.sqrt(1 - cga * cga);
  return [as, bs * cga, cs * cbe, 0, bs * sga, -cs * sbe * cal, 0, 0, 1 / c];
}

/** Direct cell from a UB in lambda units (divide by the wavelength first, or pass lambda). */
export function cellFromUB(ub, lambda = 1) {
  const u = ub.map((v) => v / lambda);
  const g = inv3(mul(transpose(u), u));
  const [a, b, c] = [0, 4, 8].map((i) => Math.sqrt(g[i]));
  const ang = (x, p, q) => Math.acos(x / (p * q)) / DEG;
  return [a, b, c, ang(g[5], b, c), ang(g[2], a, c), ang(g[1], a, b)];
}

/** Orthogonal U with UB ~ U B lambda (polar decomposition by Newton iteration). */
export function polarU(ub, cell, lambda) {
  const B = bMatrix(...cell).map((v) => v * lambda);
  let U = mul(ub, inv3(B));
  for (let k = 0; k < 30; k++) U = U.map((v, i) => 0.5 * (v + inv3(transpose(U))[i]));
  return U;
}

/** Cell parameters refined for each crystal system, and how to expand them. */
export const CELL_SYSTEMS = {
  cubic: { names: ['a'], cell: ([a]) => [a, a, a, 90, 90, 90] },
  hexagonal: { names: ['a', 'c'], cell: ([a, c]) => [a, a, c, 90, 90, 120] },
  tetragonal: { names: ['a', 'c'], cell: ([a, c]) => [a, a, c, 90, 90, 90] },
  orthorhombic: { names: ['a', 'b', 'c'], cell: ([a, b, c]) => [a, b, c, 90, 90, 90] },
  monoclinic: { names: ['a', 'b', 'c', 'beta'], cell: ([a, b, c, be]) => [a, b, c, 90, be, 90] },
  triclinic: { names: ['a', 'b', 'c', 'alpha', 'beta', 'gamma'], cell: (v) => v.slice(0, 6) },
};

/** Crystal system from a CrysAlis Laue class ("6/m", "-3m", "4/mmm", ...) and a starting cell. */
export function crystalSystem(laue, cell) {
  const near = (x, y) => Math.abs(x - y) < 1;
  const [a, b, c, al, be, ga] = cell;
  const want = /6|-3|3/.test(laue ?? '') ? 'hexagonal' : /m-3|^23/.test(laue ?? '') ? 'cubic'
    : /4/.test(laue ?? '') ? 'tetragonal' : /mmm/.test(laue ?? '') ? 'orthorhombic' : /2\/m/.test(laue ?? '') ? 'monoclinic' : 'triclinic';
  const fits = {
    hexagonal: near(al, 90) && near(be, 90) && near(ga, 120) && Math.abs(a - b) / a < 0.01,
    cubic: near(al, 90) && near(be, 90) && near(ga, 90) && Math.abs(a - b) / a < 0.01 && Math.abs(a - c) / a < 0.01,
    tetragonal: near(al, 90) && near(be, 90) && near(ga, 90) && Math.abs(a - b) / a < 0.01,
    orthorhombic: near(al, 90) && near(be, 90) && near(ga, 90),
    monoclinic: near(al, 90) && near(ga, 90),
    triclinic: true,
  };
  return fits[want] ? want : 'triclinic';
}

// ---- Instrument model -----------------------------------------------------------

/** Instrument model from a frame header (parseRodHeader), with the software zero corrections. */
export function headerGeometry(h) {
  return {
    wavelength: h.wavelengths.alpha12, ox: h.origin[0], oy: h.origin[1], distance: h.distance,
    pixelMM: h.pixelMM[0] * h.binning[0], d1: h.detRot[0], d2: h.detRot[1],
    theta: h.start[1], thetaOffset: h.zeroCorr[1], alpha: h.alpha,
    omegaOffset: h.zeroCorr[0], kappaOffset: h.zeroCorr[2], phiOffset: h.zeroCorr[3],
    tx: 0, ty: 0, tz: 0,
  };
}

/** Precomputed matrices of an instrument model `g` (see headerGeometry for the fields). */
export function prepare(g) {
  const Rdet = mul(mul(rot([0, -1, 0], g.theta + g.thetaOffset), rot([-1, 0, 0], g.d2)), rot([0, 0, 1], g.d1));
  const a = g.alpha * DEG;
  return {
    g, Rdet, normal: mulv(Rdet, [0, 0, -1]),
    eOmega: [0, -1, 0], eKappa: [0, -Math.cos(a), Math.sin(a)], ePhi: [0, -1, 0],
  };
}

/** Lab position (mm) of the centre of pixel (i, j). */
export function pixelToLab(G, i, j) {
  const { g, Rdet } = G;
  return mulv(Rdet, [(i - g.ox) * g.pixelMM, (j - g.oy) * g.pixelMM, -g.distance]);
}

/** Goniometer rotation for header angles (degrees). */
export function gonio(G, omega, kappa, phi) {
  const { g } = G;
  return mul(mul(rot(G.eOmega, omega + g.omegaOffset), rot(G.eKappa, kappa + g.kappaOffset)), rot(G.ePhi, phi + (g.phiOffset ?? 0)));
}

/** Crystal position in the lab (mm): the offset (tx, ty, tz) rides on the phi axis. */
export const samplePosition = (G, R) => mulv(R, [G.g.tx, G.g.ty, G.g.tz]);

/** Pixel (i, j) hit by a ray from `origin` (lab mm) along `s`; null if it points away. */
export function labToPixel(G, s, origin = [0, 0, 0]) {
  const { g, Rdet, normal } = G;
  const sn = dot(s, normal);
  if (!(sn > 0)) return null;
  const t = (g.distance - dot(origin, normal)) / sn;
  const P = [origin[0] + s[0] * t, origin[1] + s[1] * t, origin[2] + s[2] * t];
  const loc = mulv(transpose(Rdet), P);
  return [loc[0] / g.pixelMM + g.ox, loc[1] / g.pixelMM + g.oy];
}

/**
 * Predict where reciprocal vector xC (lambda units, phi-axis frame) diffracts during a
 * scan about `axis` (0 = omega, 3 = phi) at fixed other header angles, choosing the
 * crossing nearest `near` (header degrees). Returns { angle, i, j } or null.
 */
export function predict(G, xC, angles, axis, near) {
  const [omega, , kappa, phi] = angles;
  const s0 = [0, 0, -1];
  // x_lab = A R(e, a) v: e is the scan axis in its own frame, A the fixed part in front.
  let A, v, e, base;
  if (axis === 0) {
    A = [1, 0, 0, 0, 1, 0, 0, 0, 1];
    v = mulv(mul(rot(G.eKappa, kappa + G.g.kappaOffset), rot(G.ePhi, phi + (G.g.phiOffset ?? 0))), xC);
    e = G.eOmega;
    base = G.g.omegaOffset;
  } else {
    A = mul(rot(G.eOmega, omega + G.g.omegaOffset), rot(G.eKappa, kappa + G.g.kappaOffset));
    v = xC;
    e = G.ePhi;
    base = G.g.phiOffset ?? 0;
  }
  const t = mulv(transpose(A), s0);
  const vpar = dot(v, e), vperp = [v[0] - vpar * e[0], v[1] - vpar * e[1], v[2] - vpar * e[2]];
  const P = dot(vperp, t), Q = dot(cross(e, vperp), t);
  const C = -0.5 * dot(v, v) - vpar * dot(e, t);
  const R = Math.hypot(P, Q);
  if (!(Math.abs(C) <= R)) return null;
  const b0 = Math.atan2(Q, P), d = Math.acos(C / R);
  let best = null;
  for (const sg of [1, -1]) {
    const ang = (b0 + sg * d) / DEG - base;
    const diff = ((ang - near) % 360 + 540) % 360 - 180;
    if (!best || Math.abs(diff) < Math.abs(best)) best = diff;
  }
  const angle = near + best;
  const R1 = axis === 0 ? rot(e, angle + base) : mul(A, rot(e, angle + base));
  const x = mulv(R1, v);
  const s = [x[0] + s0[0], x[1] + s0[1], x[2] + s0[2]];
  const full = axis === 0 ? gonio(G, angle, kappa, phi) : gonio(G, omega, kappa, angle);
  const ij = labToPixel(G, s, samplePosition(G, full));
  return ij ? { angle, i: ij[0], j: ij[1] } : null;
}

/** Unit diffracted direction and scattering vector of pixel (i, j) at goniometer rotation R. */
export function pixelX(G, i, j, R) {
  const P = pixelToLab(G, i, j), c = samplePosition(G, R);
  const D = [P[0] - c[0], P[1] - c[1], P[2] - c[2]], n = norm(D);
  return [D[0] / n, D[1] / n, D[2] / n + 1];
}

// ---- Least squares ---------------------------------------------------------------

/**
 * Levenberg-Marquardt with a numerical Jacobian and a robust soft-L1 loss (iteratively
 * reweighted, scale `fScale` in residual units). `f(p)` returns a Float64Array of
 * residuals; `steps` are the finite-difference steps. Returns { p, cost, sigma, iterations }.
 */
export function levenbergMarquardt(f, p0, { steps, maxIter = 60, fScale = 1.5, onIter } = {}) {
  let p = Float64Array.from(p0);
  const n = p.length;
  let r = f(p), lambda = 1e-3;
  const weights = (res) => Float64Array.from(res, (x) => 1 / Math.sqrt(Math.sqrt(1 + (x / fScale) ** 2)));
  const costOf = (res) => res.reduce((s, x) => s + 2 * fScale * fScale * (Math.sqrt(1 + (x / fScale) ** 2) - 1), 0);
  let cost = costOf(r), iter = 0, J = null;
  for (; iter < maxIter; iter++) {
    const m = r.length, w = weights(r);
    J = Array.from({ length: n }, () => new Float64Array(m));
    for (let k = 0; k < n; k++) {
      const q = Float64Array.from(p);
      q[k] += steps[k];
      const rk = f(q);
      for (let i = 0; i < m; i++) J[k][i] = (rk[i] - r[i]) / steps[k];
    }
    const A = new Float64Array(n * n), g = new Float64Array(n);
    for (let a = 0; a < n; a++) {
      for (let i = 0; i < m; i++) g[a] += J[a][i] * w[i] * w[i] * r[i];
      for (let b = a; b < n; b++) {
        let s = 0;
        for (let i = 0; i < m; i++) s += J[a][i] * w[i] * w[i] * J[b][i];
        A[a * n + b] = A[b * n + a] = s;
      }
    }
    let improved = false;
    for (let tries = 0; tries < 12; tries++) {
      const Ad = Float64Array.from(A);
      for (let a = 0; a < n; a++) Ad[a * n + a] += lambda * (A[a * n + a] || 1e-12);
      const dp = solve(Ad, g.map((x) => -x), n);
      if (!dp) { lambda *= 10; continue; }
      const q = p.map((x, k) => x + dp[k]);
      const rq = f(q), cq = costOf(rq);
      if (cq < cost) {
        const rel = (cost - cq) / Math.max(cost, 1e-300);
        p = q; r = rq; cost = cq; lambda = Math.max(lambda / 10, 1e-9); improved = true;
        onIter?.(iter, cost);
        if (rel < 1e-10) iter = maxIter;
        break;
      }
      lambda *= 10;
    }
    if (!improved) break;
  }
  // standard deviations from the (weighted) normal matrix at the solution
  const m = r.length, w = weights(r);
  const A = new Float64Array(n * n);
  for (let a = 0; a < n; a++) {
    for (let b = a; b < n; b++) {
      let s = 0;
      for (let i = 0; i < m; i++) s += J[a][i] * w[i] * w[i] * J[b][i];
      A[a * n + b] = A[b * n + a] = s;
    }
  }
  const s2 = r.reduce((s, x, i) => s + (w[i] * x) ** 2, 0) / Math.max(1, m - n);
  const sigma = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const e = new Float64Array(n);
    e[k] = 1;
    const col = solve(Float64Array.from(A), e, n);
    sigma[k] = col ? Math.sqrt(Math.abs(col[k]) * s2) : NaN;
  }
  return { p, cost, residuals: r, sigma, iterations: iter };
}

/** Solve A x = b (n x n, row-major Float64Array) by Gaussian elimination with pivoting. */
function solve(A, b, n) {
  const x = Float64Array.from(b);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r * n + c]) > Math.abs(A[piv * n + c])) piv = r;
    if (!(Math.abs(A[piv * n + c]) > 1e-300)) return null;
    if (piv !== c) {
      for (let k = 0; k < n; k++) [A[c * n + k], A[piv * n + k]] = [A[piv * n + k], A[c * n + k]];
      [x[c], x[piv]] = [x[piv], x[c]];
    }
    for (let r = c + 1; r < n; r++) {
      const f = A[r * n + c] / A[c * n + c];
      if (!f) continue;
      for (let k = c; k < n; k++) A[r * n + k] -= f * A[c * n + k];
      x[r] -= f * x[c];
    }
  }
  for (let r = n - 1; r >= 0; r--) {
    let s = x[r];
    for (let k = r + 1; k < n; k++) s -= A[r * n + k] * x[k];
    x[r] = s / A[r * n + r];
  }
  return x;
}
