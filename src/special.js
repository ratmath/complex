// The gamma and zeta functions on the complex plane, as proven enclosures.
//
// Each is one formula with a bounded remainder, written once and run either
// on rectangles (a value at a point or over an enclosure) or on Taylor models
// (the function on a disc, ready for contour integrals, the argument
// principle and zero finding).
//
//   gamma: the reflection formula left of Re z = 1/2, then the recurrence
//          to push z right, then Stirling's series, whose remainder is at
//          most sec^2n(arg z / 2) times the first term left out.
//   zeta:  the functional equation left of Re s = 0, otherwise the
//          Euler–Maclaurin formula, whose remainder after the term in B_2m
//          is at most |B_2m|/(2m)! · |s(s+1)…(s+2m−1)| · N^(1−σ−2m)/(σ+2m−1).
//
// Both are computed through their entire companions 1/Γ(z) and (s − 1)ζ(s),
// so that nothing is divided by until the very end.

import { Rational, RationalInterval } from "@ratmath/core";
import { enclose, constants, bernoulli } from "@ratmath/analytic";
import { Complex, Q, ZERO, ONE, TWO, refuse } from "./complex.js";
import { Rect, Disc, I, iv, sqrtUp, sqrtDown } from "./enclosure.js";
import { exp, log, sin } from "./functions.js";
import { CTM, _registerModel } from "./taylor.js";
import { zeros } from "./contour.js";

const HALF = new Rational(1n, 2n);
const R = (n, d = 1n) => new Rational(BigInt(n), BigInt(d));
const isModel = (x) => x instanceof CTM;
/** A rectangle that holds every value of x. */
const hullOf = (x) => (isModel(x) ? Rect.from(x.range()) : x);
/** x with an error of size at most b allowed for. */
const widen = (x, b) => (isModel(x) ? x._more(b) : x.add(new Rect(iv(b.negate(), b), iv(b.negate(), b))));
const ex = (x, bits) => (isModel(x) ? exp(x) : exp(x.round(bits + 16), { bits }));
const ln = (x, bits) => (isModel(x) ? log(x) : log(x.round(bits + 16), { bits }).rect);
const realRect = (i) => new Rect(iv(i.low, i.high), iv(ZERO, ZERO));
const piRect = (bits) => realRect(constants.pi(bits));
const lnCache = new Map();
/** ln n for a whole number n, as a rectangle on the real axis. */
const lnOf = (n, bits) => {
  const key = `${n}/${bits}`;
  if (!lnCache.has(key)) {
    if (lnCache.size > 20000) lnCache.clear();
    lnCache.set(key, realRect(enclose.log(R(n), bits)));
  }
  return lnCache.get(key);
};
const oneMinus = (x) => x.neg().add(Complex.ONE);
/** Rectangles are exact, so their fractions grow with every product: round outward to the working precision. */
const trim = (x, bits) => (isModel(x) ? x : x.round(bits + 16));
const plus = (x, q) => x.add(new Complex(q, ZERO));
const eps = (bits) => new Rational(1n, 1n << BigInt(bits + 8));
const absRange = (H) => {
  const n = H.normRange();
  return { lo: sqrtDown(n.lo), hi: sqrtUp(n.hi) };
};

// ---- gamma ----------------------------------------------------------------------------------

/** log Γ(y) for Re y well to the right, by Stirling's series with its remainder; null when y is not far enough right. */
function stirling(y, bits) {
  const H = hullOf(y);
  const { lo: m, hi: M } = absRange(H);
  const reLo = H.re.lo;
  if (reLo.numerator <= 0n) return null;
  // sec²(θ/2) = 2|y|/(|y| + Re y)
  const sec2 = TWO.multiply(M).divide(m.add(reLo));
  const target = eps(bits);
  let n = 1;
  let bound = null;
  let power = m; // m^(2n−1)
  let sec = sec2; // sec2^n
  for (; n <= 400; n++) {
    bound = bernoulli(2 * n).abs().divide(R(2 * n * (2 * n - 1))).multiply(sec).divide(power);
    if (bound.compareTo(target) < 0) break;
    power = power.multiply(m).multiply(m);
    sec = sec.multiply(sec2);
  }
  if (n > 400) return null;
  const w = trim(y.inv(), bits);
  const w2 = trim(w.mul(w), bits);
  let S = null;
  for (let k = n - 1; k >= 1; k--) {
    const c = bernoulli(2 * k).divide(R(2 * k * (2 * k - 1)));
    S = S === null ? plus(w2.scale(ZERO), c) : plus(trim(S.mul(w2), bits), c);
  }
  const series = S === null ? null : trim(S.mul(w), bits);
  const pi = constants.pi(bits + 8);
  const half2pi = realRect(enclose.log(new RationalInterval(pi.low.multiply(TWO), pi.high.multiply(TWO)), bits + 8)).scale(HALF);
  let L = y.sub(new Complex(HALF, ZERO)).mul(ln(y, bits + 8)).sub(y).add(half2pi);
  if (series) L = L.add(series);
  return widen(L, bound);
}

/**
 * The shift that takes x far enough right for Stirling's series, and log Γ
 * there: { N, L } with L = log Γ(x + N). The enclosure must lie in Re > 0.
 */
function shifted(x, bits, name) {
  const H = hullOf(x);
  if (H.re.lo.numerator <= 0n) refuse("precision", `${name}: the enclosure is too wide (it straddles Re z = 1/2 and reaches Re z ≤ 0); use a smaller one`);
  for (let M = Math.ceil(bits / 6) + 3; M < 4 * bits + 64; M *= 2) {
    const shift = R(M).subtract(H.re.lo);
    const N = shift.numerator > 0n ? Number(shift.numerator / shift.denominator) + 1 : 0;
    const L = stirling(plus(x, R(N)), bits);
    if (L) return { N, L };
  }
  return refuse("precision", `${name}: the working precision could not be reached; use a smaller enclosure`);
}
/** log Γ(x) for an enclosure in Re > 0: log Γ(x + N) − Σ log(x + k). */
function lgammaRight(x, bits, name = "gamma") {
  let { N, L } = shifted(x, bits, name);
  for (let k = 0; k < N; k++) L = L.sub(ln(plus(x, R(k)), bits + 8));
  return L;
}
/** 1/Γ(x) for an enclosure in Re > 0: exp(−log Γ(x + N)) · x(x + 1)…(x + N − 1). */
function rgammaRight(x, bits) {
  const { N, L } = shifted(x, bits, "gamma");
  if (!isModel(x)) {
    let out = ex(L.neg(), bits + 8);
    for (let k = 0; k < N; k++) out = trim(out.mul(plus(x, R(k))), bits);
    return out;
  }
  // a model works to a fixed number of places, and a tiny exponential times a long product of large factors
  // would lose them to the difference in size: each factor is divided by a whole number near its size, and the
  // logarithms of those numbers go into the exponent
  const c = hullOf(x).center;
  const sizes = [];
  let E = L.neg();
  for (let k = 0; k < N; k++) {
    const [u, v] = plus(c, R(k)).toNumbers();
    const size = Math.max(1, Math.round(Math.hypot(u, v)));
    sizes.push(size);
    if (size > 1) E = E.add(lnOf(size, bits + 8));
  }
  let out = ex(E, bits + 8);
  for (let k = 0; k < N; k++) out = out.mul(plus(x, R(k)).scale(R(1, sizes[k])));
  return out;
}
const sinPi = (x, bits) => sin(trim(x.mul(piRect(bits + 8)), bits), isModel(x) ? undefined : { bits: bits + 8 });
/** 1/Γ(x), entire. x is a Rect or a model. */
function rgammaCore(x, bits) {
  const re = hullOf(x).center.re;
  if (re.compareTo(HALF) >= 0) return rgammaRight(x, bits);
  if (re.compareTo(R(-64)) > 0) {
    // 1/Γ(x) = x(x + 1)…(x + K − 1) / Γ(x + K): products only
    const gap = HALF.subtract(re);
    const K = Number(gap.numerator / gap.denominator) + 1;
    let out = rgammaRight(plus(x, R(K)), bits);
    for (let k = 0; k < K; k++) out = trim(out.mul(plus(x, R(k))), bits);
    return out;
  }
  // far to the left, the reflection formula: 1/Γ(x) = sin(πx)·Γ(1 − x)/π
  return trim(sinPi(x, bits).mul(gammaCore(oneMinus(x), bits)), bits).mul(piRect(bits + 8).inv());
}
/** Γ(x). */
function gammaCore(x, bits) {
  if (hullOf(x).center.re.compareTo(HALF) >= 0) {
    return rgammaRight(x, bits).inv();
  }
  return rgammaCore(x, bits).inv();
}
const holdsPole = (H) => {
  // a whole number ≤ 0 inside the rectangle
  if (!I.hasZero(H.im) || H.re.lo.numerator > 0n) return false;
  const hi = H.re.hi.numerator > 0n ? ZERO : H.re.hi;
  const k = hi.numerator / hi.denominator; // truncates towards 0: the largest whole number ≤ hi when hi ≤ 0 is this or one less
  const top = R(k).compareTo(hi) <= 0 ? R(k) : R(k - 1n);
  return top.compareTo(H.re.lo) >= 0;
};

const prepare = (z, opts) => {
  const bits = typeof opts === "number" ? opts : opts?.bits ?? 64;
  if (isModel(z)) return { x: z, bits: z.bits - 16, done: (v) => v };
  const x = Rect.from(z);
  if (x.isPoint) return { x, bits: bits + 24, done: (v) => v.round(bits + 4) };
  // over an enclosure, a Taylor model on the disc around it keeps the result about as wide as the function
  // really varies; plain interval arithmetic would add up the width of every step
  const D = x.toDisc();
  return { x: CTM.z(D, { order: 8, bits: bits + 56 }), bits: bits + 24, done: (v) => Rect.from(v.range()).round(bits + 16) };
};

/** 1/Γ(z): an entire function, with zeros at 0, −1, −2, … */
export function rgamma(z, opts) {
  const { x, bits, done } = prepare(z, opts);
  return done(rgammaCore(x, bits));
}
/**
 * Γ(z) for complex z: a Complex, a Rect or a Disc gives a Rect certain to
 * contain every value; a Taylor model gives the model of Γ of it. Refuses at
 * the poles 0, −1, −2, …
 */
export function gamma(z, opts) {
  const { x, bits, done } = prepare(z, opts);
  if (holdsPole(hullOf(x))) refuse("pole", "gamma: the input may be one of the poles 0, −1, −2, …");
  return done(gammaCore(x, bits));
}
/**
 * log Γ(z) for Re z > 0, the branch that is real on the positive real axis
 * (it is not always the principal logarithm of Γ(z)).
 */
export function logGamma(z, opts) {
  const { x, bits, done } = prepare(z, opts);
  if (hullOf(x).re.lo.numerator <= 0n) refuse("bad-argument", "logGamma is for Re z > 0; elsewhere use gamma");
  return done(lgammaRight(x, bits, "logGamma"));
}

// ---- zeta -----------------------------------------------------------------------------------

/** least[n]: the smallest prime factor of n, for n up to N. */
function leastFactors(N) {
  const least = new Int32Array(N + 1);
  for (let p = 2; p <= N; p++) {
    if (least[p]) continue;
    for (let q = p; q <= N; q += p) if (!least[q]) least[q] = p;
  }
  return least;
}
/** (s − 1)·ζ(s), entire. */
function zeta1Core(s, bits) {
  const H = hullOf(s);
  if (H.re.hi.numerator < 0n) {
    // ζ(s) = 2^s π^(s−1) sin(πs/2) Γ(1 − s) ζ(1 − s), and ζ(1 − s) = [(1 − s) − 1]ζ(1 − s) / (−s)
    const t = oneMinus(s);
    const pi = piRect(bits + 8);
    const lnPi = realRect(enclose.log(new RationalInterval(constants.pi(bits + 8).low, constants.pi(bits + 8).high), bits + 8));
    const z = ex(s.mul(lnOf(2, bits + 8)), bits + 8)
      .mul(ex(s.sub(Complex.ONE).mul(lnPi), bits + 8))
      .mul(sin(s.mul(pi.scale(HALF)), isModel(s) ? undefined : { bits: bits + 8 }))
      .mul(gammaCore(t, bits))
      .mul(zeta1Core(t, bits))
      .mul(s.neg().inv());
    return z.mul(s.sub(Complex.ONE));
  }
  const { hi: A } = absRange(H);
  const sigma = H.re.lo;
  // the whole number at or below σ, for powers of N
  let fl = sigma.numerator / sigma.denominator;
  if (R(fl).compareTo(sigma) > 0) fl -= 1n;
  const target = eps(bits);
  const s1Abs = absRange(H.sub(Complex.ONE)).hi;
  for (let N = Math.max(8, Math.ceil(A.toNumber() / 2 + bits / 4)); N < 1 << 22; N *= 2) {
    // the remainder after the term in B_2m: |B_2m|/(2m)! · Π_{j<2m}(|s| + j) · N^(1−σ−2m)/(σ + 2m − 1)
    let m = 1;
    let bound = null;
    let rising = A.multiply(A.add(ONE)); // Π_{j<2m}(A + j)
    let fact = 2n; // (2m)!
    for (; m <= 150; m++) {
      const denom = sigma.add(R(2 * m - 1));
      if (denom.numerator > 0n) {
        const e = 1n - fl - BigInt(2 * m); // N^(1 − ⌊σ⌋ − 2m) ≥ N^(1 − σ − 2m)
        const power = e >= 0n ? R(BigInt(N) ** e) : R(1n, BigInt(N) ** -e);
        bound = bernoulli(2 * m).abs().divide(R(fact)).multiply(rising).multiply(power).divide(denom).multiply(s1Abs);
        if (bound.compareTo(target) < 0) break;
      }
      rising = rising.multiply(A.add(R(2 * m))).multiply(A.add(R(2 * m + 1)));
      fact *= BigInt(2 * m + 1) * BigInt(2 * m + 2);
    }
    if (m > 150) continue;
    // Σ_{n<N} n^−s
    // n^−s is an exponential for a prime n and a product of earlier ones otherwise
    const least = leastFactors(N);
    const pw = [null, plus(s.scale(ZERO), ONE)];
    let sum = pw[1];
    for (let n = 2; n <= N; n++) {
      pw[n] = least[n] === n ? ex(trim(s.mul(lnOf(n, bits + 8)), bits).neg(), bits + 8) : trim(pw[least[n]].mul(pw[n / least[n]]), bits);
      if (n < N) sum = sum.add(pw[n]);
    }
    const Ns = pw[N]; // N^−s
    // Σ_{k≤m} U_k with U_k = B_2k/(2k)! · s(s+1)…(s+2k−2)/N^(2k−1), each made from the one before so that
    // no factor is ever far from the size of the term itself
    const N2 = R(BigInt(N) * BigInt(N));
    let U = s.scale(bernoulli(2).divide(R(2 * N)));
    let corr = plus(s.scale(ZERO), HALF).add(U);
    for (let k = 1; k < m; k++) {
      const ratio = bernoulli(2 * k + 2).divide(bernoulli(2 * k)).divide(R((2 * k + 1) * (2 * k + 2))).divide(N2);
      U = trim(trim(U.mul(plus(s, R(2 * k - 1))), bits).mul(plus(s, R(2 * k))), bits).scale(ratio);
      corr = corr.add(U);
    }
    // (s − 1)ζ(s) = (s − 1)[Σ + N^−s(1/2 + corrections)] + N^(1−s)
    const body = trim(trim(sum.add(Ns.mul(corr)), bits).mul(s.sub(Complex.ONE)), bits).add(Ns.scale(R(N)));
    return widen(body, bound);
  }
  return refuse("precision", "zeta: the working precision could not be reached; use a smaller enclosure");
}

/** (s − 1)·ζ(s): an entire function. */
export function zetaTimesSMinus1(s, opts) {
  const { x, bits, done } = prepare(s, opts);
  return done(zeta1Core(x, bits));
}
/**
 * The Riemann zeta function for complex s: a Complex, a Rect or a Disc gives
 * a Rect certain to contain every value; a Taylor model gives the model of ζ
 * of it. Refuses at the pole s = 1.
 */
export function zeta(s, opts) {
  const { x, bits, done } = prepare(s, opts);
  const H = hullOf(x);
  if (H.contains(Complex.ONE)) refuse("pole", "zeta: the input may be the pole s = 1");
  return done(zeta1Core(x, bits).mul(x.sub(Complex.ONE).inv()));
}

_registerModel("gamma", (f) => gamma(f));
_registerModel("rgamma", (f) => rgamma(f));
_registerModel("lgamma", (f) => logGamma(f));
_registerModel("zeta", (f) => zeta(f));
_registerModel("zeta1", (f) => zetaTimesSMinus1(f));

/**
 * The zeros of ζ with imaginary part between t1 and t2 (the non-trivial ones:
 * they all lie in 0 < Re s < 1).
 *
 *   zetaZeros(10, 30)
 *   → { zeros: [{ t, enclosure, onCriticalLine }], count, complete, allOnCriticalLine }
 *
 * Each zero is isolated and proved simple by `zeros`, and `t` is the interval
 * holding its imaginary part. Then the Riemann hypothesis is checked for it:
 * a small disc centred on the critical line is shown to hold exactly one
 * zero. Zeros come in pairs s and 1 − s̄, mirror images in that line, and the
 * disc is its own mirror image; a single zero in it is therefore its own
 * mirror image, which puts its real part at exactly 1/2.
 */
export function zetaZeros(t1, t2, opts = {}) {
  const [a, b] = [Q(t1), Q(t2)];
  if (a.compareTo(b) >= 0 || a.numerator <= 0n) throw new RangeError("give 0 < t1 < t2");
  const bits = opts.bits ?? 64;
  const found = zeros((z) => zeta(z), new Rect(iv(ZERO, ONE), iv(a, b)), { order: opts.order ?? 10, bits, tolerance: opts.tolerance, maxCells: opts.maxCells });
  const out = [];
  for (const z of found.zeros) {
    let onCriticalLine = false;
    if (z.simple) {
      const c = z.enclosure.center;
      const off = c.re.subtract(HALF).abs().add(z.enclosure.size);
      const base = new Rational(1n, 1n << 20n);
      const rho = off.multiply(R(4)).compareTo(base) > 0 ? off.multiply(R(4)) : base;
      try {
        const D = new Disc(new Complex(HALF, Rect.from(new Complex(ZERO, c.im)).round(40).center.im), rho);
        const m = zeta(CTM.z(D, { order: 6, bits: 96 + 2 * 20 }));
        onCriticalLine = m.zeroCount() === 1 && D.contains(z.enclosure);
      } catch (e) {
        if (!e.code) throw e;
      }
    }
    out.push({ t: z.enclosure.im, enclosure: z.enclosure, multiplicity: z.multiplicity, onCriticalLine, boundary: z.boundary });
  }
  out.sort((p, q) => p.t.lo.compareTo(q.t.lo));
  return { zeros: out, count: found.count, complete: found.complete, allOnCriticalLine: out.every((z) => z.onCriticalLine), unresolved: found.unresolved };
}
