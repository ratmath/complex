// Complex Taylor models: a polynomial plus a proven bound on what it leaves out.
//
// A model of f on the closed disc |z − c| ≤ r is a polynomial P of degree at
// most n and a number ε with
//
//     |f(z) − P(z − c)| ≤ ε    for every z in the disc,
//
// and f is analytic there. Every operation keeps that statement true: sums and
// products push the terms they drop into ε, and exp, log, roots and quotients
// are built from their series with the tail bounded. When a hypothesis fails
// (a divisor that may vanish on the disc, a logarithm whose argument may be 0,
// a series that converges too slowly to bound) the operation refuses instead
// of returning something unproved, and the caller takes a smaller disc.
//
// Inside, the polynomial is in the scaled variable w = (z − c)/r, |w| ≤ 1, so
// that the size of P on the disc is at most the sum of the sizes of its
// coefficients. Coefficients are fixed-point: pairs of BigInts counting units
// of 2^−bits. Rounding is always towards −∞ and every rounding adds to ε,
// which is kept as a whole number of the same units.

import { Rational } from "@ratmath/core";
import { Complex, Q, ZERO, ONE, isqrt, refuse } from "./complex.js";
import { Rect, Disc, iv } from "./enclosure.js";
import { exp as expPoint, log as logPoint, root as rootPoint, _setModelHook } from "./functions.js";
import { constants } from "@ratmath/analytic";

const floorDiv = (a, b) => {
  const q = a / b;
  return a % b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
};
const ceilDiv = (a, b) => -floorDiv(-a, b);
/** An upper bound, in units, for the modulus of a fixed-point pair. */
const up = ([x, y]) => (x === 0n && y === 0n ? 0n : y === 0n ? (x < 0n ? -x : x) : x === 0n ? (y < 0n ? -y : y) : isqrt(x * x + y * y) + 1n);
const down = ([x, y]) => isqrt(x * x + y * y);
const ZERO2 = [0n, 0n];

export class CTM {
  /** Use taylorModel(f, disc) or CTM.z(disc). */
  constructor(center, radius, a, e, order, P) {
    this.center = center;
    this.radius = radius;
    this.a = a;
    this.e = e;
    this.order = order;
    this.P = P;
    Object.freeze(this);
  }
  static _opts(opts = {}) {
    const order = opts.order ?? 12;
    const bits = opts.bits ?? 96;
    if (!Number.isInteger(order) || order < 1 || order > 200) throw new RangeError("order must be a whole number from 1 to 200");
    if (!Number.isInteger(bits) || bits < 16 || bits > 4000) throw new RangeError("bits must be a whole number from 16 to 4000");
    return { order, P: BigInt(bits) };
  }
  static _disc(d) {
    const D = d instanceof Disc ? d : d instanceof Rect ? d.toDisc() : Disc.from(d);
    if (D.radius.numerator <= 0n) throw new RangeError("a Taylor model needs a disc of positive radius");
    return D;
  }
  /** The model of the function z itself on a disc. Options: order (default 12), bits (default 96). */
  static z(d, opts) {
    const D = CTM._disc(d);
    const { order, P } = CTM._opts(opts);
    const one = 1n << P;
    const c = D.center;
    const a0 = [floorDiv(c.re.numerator * one, c.re.denominator), floorDiv(c.im.numerator * one, c.im.denominator)];
    const a1 = [floorDiv(D.radius.numerator * one, D.radius.denominator), 0n];
    return new CTM(c, D.radius, [a0, a1], 3n, order, P);
  }
  /** A constant as a model on a disc: an exact number, or a Rect or Disc enclosing one. */
  static constant(d, value, opts) {
    const D = CTM._disc(d);
    const { order, P } = CTM._opts(opts);
    return liftTo({ center: D.center, radius: D.radius, order, P }, value);
  }
  get bits() {
    return Number(this.P);
  }
  get disc() {
    return new Disc(this.center, this.radius);
  }
  _like(a, e, order = this.order) {
    while (a.length > 1 && a[a.length - 1][0] === 0n && a[a.length - 1][1] === 0n) a.pop();
    return new CTM(this.center, this.radius, a, e, order, this.P);
  }
  _same(o) {
    if (!o.center.equals(this.center) || !o.radius.equals(this.radius) || o.P !== this.P) {
      refuse("bad-argument", "the two Taylor models are on different discs or at different precisions");
    }
  }
  _lift(x) {
    if (x instanceof CTM) {
      this._same(x);
      return x;
    }
    return liftTo(this, x);
  }
  /** Σ over k ≥ from of the bounds of the coefficients, in units. */
  _sum(from = 0) {
    let s = 0n;
    for (let k = from; k < this.a.length; k++) s += up(this.a[k]);
    return s;
  }
  /** An upper bound for |f − f(c)-part|: the size of everything but the constant term, with ε, as a Rational. */
  _spread() {
    return new Rational(this._sum(1) + this.e, 1n << this.P);
  }
  _q(n) {
    return new Rational(n, 1n << this.P);
  }
  _c(pair) {
    return new Complex(this._q(pair[0]), this._q(pair[1]));
  }

  // ---- arithmetic ----

  add(o) {
    o = this._lift(o);
    const n = Math.max(this.a.length, o.a.length);
    const a = [];
    for (let k = 0; k < n; k++) {
      const [x, y] = [this.a[k] ?? ZERO2, o.a[k] ?? ZERO2];
      a.push([x[0] + y[0], x[1] + y[1]]);
    }
    return this._like(a, this.e + o.e, Math.min(this.order, o.order));
  }
  neg() {
    // −floor(x) is within one unit of floor(−x): the same model, exactly negated
    return this._like(this.a.map(([x, y]) => [-x, -y]), this.e);
  }
  sub(o) {
    return this.add(this._lift(o).neg());
  }
  mul(o) {
    o = this._lift(o);
    const order = Math.min(this.order, o.order);
    const P = this.P;
    const [A, B] = [this.a, o.a];
    const len = A.length + B.length - 1;
    const re = new Array(len).fill(0n);
    const im = new Array(len).fill(0n);
    for (let i = 0; i < A.length; i++) {
      const [ar, ai] = A[i];
      if (ar === 0n && ai === 0n) continue;
      for (let j = 0; j < B.length; j++) {
        const [br, bi] = B[j];
        re[i + j] += ar * br - ai * bi;
        im[i + j] += ar * bi + ai * br;
      }
    }
    const a = [];
    // each kept coefficient is rounded down: at most one unit in each part, under 2 in modulus
    let e = 0n;
    for (let k = 0; k < len; k++) {
      const c = [re[k] >> P, im[k] >> P];
      if (k <= order) {
        a.push(c);
        e += 2n;
      } else e += up(c) + 2n;
    }
    const [SA, SB] = [this._sum(), o._sum()];
    e += ceilDiv(this.e * SB + o.e * SA + this.e * o.e, 1n << P);
    return this._like(a, e, order);
  }
  /** Multiply by a real number. */
  scale(q) {
    return this.mul(Q(q));
  }
  /**
   * 1/f. Refuses when f may vanish on the disc.
   *
   * The reciprocal of the polynomial part is found as a power series (no
   * rigour needed for a candidate g), and then checked: with δ = f·g − 1
   * bounded by d < 1 on the disc, f does not vanish and
   * |1/f − g| = |g|·|δ|/|1 + δ| ≤ |g|·d/(1 − d).
   */
  inv() {
    const no = () => refuse("division-by-zero", `division by a function that may vanish on ${short(this)}`);
    const a0 = this.a[0];
    const n2 = a0[0] * a0[0] + a0[1] * a0[1];
    if (n2 === 0n) no();
    const P = this.P;
    const k0 = [(a0[0] << (2n * P)) / n2, (-a0[1] << (2n * P)) / n2];
    const g = [k0];
    for (let k = 1; k <= this.order; k++) {
      let [sr, si] = [0n, 0n];
      for (let j = 1; j <= k && j < this.a.length; j++) {
        const [x, y] = this.a[j];
        const [u, v] = g[k - j];
        sr += x * u - y * v;
        si += x * v + y * u;
      }
      sr >>= P;
      si >>= P;
      g.push([-((k0[0] * sr - k0[1] * si) >> P), -((k0[0] * si + k0[1] * sr) >> P)]);
    }
    const G = this._like(g, 0n);
    const delta = this.mul(G).sub(1);
    const d = delta._sum() + delta.e;
    const one = 1n << P;
    if (d * 8n >= one * 7n) no();
    return G._more(new Rational(G._sum() * d, (one - d) * one));
  }
  div(o) {
    return this.mul(this._lift(o).inv());
  }
  /** A whole-number power. */
  pow(n) {
    if (typeof n === "bigint") n = Number(n);
    if (!Number.isInteger(n)) return powModel(this, n);
    if (n < 0) return this.inv().pow(-n);
    let result = this._lift(1);
    let base = this;
    while (n > 0) {
      if (n & 1) result = result.mul(base);
      n >>= 1;
      if (n) base = base.mul(base);
    }
    return result;
  }
  /** The same model with its constant term removed (ε stays). */
  _tail() {
    return this._like([ZERO2, ...this.a.slice(1)], this.e);
  }
  /** The same model with a Rational added to ε. */
  _more(q) {
    return this._like(this.a.slice(), this.e + ceilDiv(q.numerator << this.P, q.denominator));
  }

  // ---- what a model says ----

  /** The proven bound ε on |f − P| over the disc, as a Rational. */
  get error() {
    return this._q(this.e);
  }
  /**
   * The Taylor polynomial's coefficients in powers of (z − c), as exact
   * Complex numbers: an approximation to f^(k)(c)/k!, for which `error`
   * together with the disc gives the guarantee.
   */
  get coefficients() {
    let rk = ONE;
    return this.a.map((p) => {
      const c = this._c(p).scale(ONE.divide(rk));
      rk = rk.multiply(this.radius);
      return c;
    });
  }
  /** A disc certain to contain f(z) for every z in the model's disc. */
  range() {
    return new Disc(this._c(this.a[0]), this._spread());
  }
  /**
   * The number of zeros of f in the disc, counted with multiplicity, when the
   * model can tell: if one coefficient outweighs all the others and ε
   * together, then by Rouché's theorem f has as many zeros in the disc as
   * that term alone, which is its index (and none on the rim). Otherwise null.
   */
  zeroCount() {
    let best = 0;
    for (let k = 1; k < this.a.length; k++) if (down(this.a[k]) > down(this.a[best])) best = k;
    let rest = this.e;
    for (let k = 0; k < this.a.length; k++) if (k !== best) rest += up(this.a[k]);
    return down(this.a[best]) > rest ? best : null;
  }
  _w(z, what) {
    z = Complex.from(z);
    const d = z.sub(this.center);
    if (d.norm().compareTo(this.radius.multiply(this.radius)) > 0) refuse("bad-argument", `${what} ${z} is outside the model's disc ${this.disc}`);
    return d.scale(ONE.divide(this.radius));
  }
  /** A rectangle certain to contain f(z), for a point z of the disc. */
  at(z) {
    const { v, err } = evalFixed(this.a, this._w(z, "the point"), this.P);
    const t = this._q(err + this.e);
    const [x, y] = [this._q(v[0]), this._q(v[1])];
    return new Rect(iv(x.subtract(t), x.add(t)), iv(y.subtract(t), y.add(t)));
  }
  /**
   * A model of f′ on the smaller disc of radius `radius` (default half). The
   * polynomial is differentiated; the part the polynomial leaves out is
   * analytic and at most ε in size, so by Cauchy's estimate its derivative is
   * at most ε/(r − radius) on the smaller disc.
   */
  derivative(radius) {
    const rho = radius === undefined ? this.radius.divide(new Rational(2n)) : Q(radius);
    if (rho.numerator <= 0n || rho.compareTo(this.radius) >= 0) throw new RangeError("the derivative lives on a strictly smaller disc: 0 < radius < the model's radius");
    const a = [];
    let e = 0n;
    let ratio = ONE.divide(this.radius); // (ρ/r)^j / r
    const step = rho.divide(this.radius);
    for (let j = 0; j + 1 < this.a.length; j++) {
      const f = ratio.multiply(new Rational(BigInt(j + 1)));
      const [x, y] = this.a[j + 1];
      a.push([floorDiv(x * f.numerator, f.denominator), floorDiv(y * f.numerator, f.denominator)]);
      e += 2n;
      ratio = ratio.multiply(step);
    }
    if (!a.length) a.push(ZERO2);
    const gap = this.radius.subtract(rho);
    e += ceilDiv(this.e * gap.denominator, gap.numerator);
    return new CTM(this.center, rho, a, e, this.order, this.P);
  }
  /**
   * A rectangle certain to contain the integral of f from a to b, two points
   * of the disc. f is analytic on the disc, so the integral does not depend
   * on the path between them.
   */
  integral(a, b) {
    [a, b] = [Complex.from(a), Complex.from(b)];
    const [wa, wb] = [this._w(a, "the endpoint"), this._w(b, "the endpoint")];
    // the antiderivative of P in w: A_{k+1} = a_k/(k + 1)
    const A = [ZERO2];
    for (let k = 0; k < this.a.length; k++) A.push([floorDiv(this.a[k][0], BigInt(k + 1)), floorDiv(this.a[k][1], BigInt(k + 1))]);
    const [Va, Vb] = [evalFixed(A, wa, this.P), evalFixed(A, wb, this.P)];
    const errs = Va.err + Vb.err + 4n * BigInt(this.a.length);
    // × r for dz = r·dw; the part left out contributes at most ε·|b − a|
    const value = new Complex(this._q(Vb.v[0] - Va.v[0]), this._q(Vb.v[1] - Va.v[1])).scale(this.radius);
    const chord = new Rational(up(fixedOf(b.sub(a), this.P)) + 2n, 1n << this.P);
    const t = this._q(errs).multiply(this.radius).add(this._q(this.e).multiply(chord));
    return new Rect(iv(value.re.subtract(t), value.re.add(t)), iv(value.im.subtract(t), value.im.add(t)));
  }
  toString(digits = 6) {
    // a part no larger than the error bound is noise from rounding, and is shown as 0
    const abs = (v) => (v < 0n ? -v : v);
    let rk = ONE;
    const cs = this.a.map(([x, y]) => {
      const c = new Complex(abs(x) > this.e ? this._q(x) : ZERO, abs(y) > this.e ? this._q(y) : ZERO).scale(ONE.divide(rk));
      rk = rk.multiply(this.radius);
      const [u, v] = c.toNumbers();
      const f = (t) => Number(t.toPrecision(digits));
      return v === 0 ? `${f(u)}` : u === 0 ? `${f(v)}i` : `(${f(u)} ${v < 0 ? "-" : "+"} ${f(Math.abs(v))}i)`;
    });
    const terms = cs.map((c, k) => (k === 0 ? c : `${c}·w${k > 1 ? `^${k}` : ""}`)).filter((_, k) => k === 0 || up(this.a[k]) > this.e);
    return `${terms.join(" + ").replace(/\+ -/g, "- ")}  ± ${Number(this.error.toNumber().toPrecision(3))}   (w = z − (${this.center}), |w| ≤ ${this.radius})`;
  }
}

/** The model's disc in a few digits, for messages. */
const short = (m) => {
  const [x, y] = m.center.toNumbers();
  const f = (v) => Number(v.toPrecision(8));
  return `the disc of radius ${f(m.radius.toNumber())} about ${f(x)}${y < 0 ? " - " : " + "}${f(Math.abs(y))}i`;
};
const fixedOf = (z, P) => [floorDiv(z.re.numerator << P, z.re.denominator), floorDiv(z.im.numerator << P, z.im.denominator)];

/** A constant (an exact number, or a Rect or Disc enclosing one) as a model shaped like `m`. */
function liftTo(m, x) {
  const P = m.P;
  let c;
  let extra = 0n;
  if (x instanceof Rect || x instanceof Disc) {
    const d = x instanceof Disc ? x : null;
    c = x.center;
    if (!x.isPoint) {
      // the half-diagonal of a rectangle, or the radius of a disc, goes into ε
      const half = d ? [floorDiv(d.radius.numerator << P, d.radius.denominator) + 1n, 0n] : fixedOf(new Complex(x.width, x.height).scale(new Rational(1n, 2n)), P).map((v) => v + 1n);
      extra = up(half);
    }
  } else c = Complex.from(x);
  const a0 = fixedOf(c, P);
  const exact = c.re.denominator === 1n && c.im.denominator === 1n;
  return new CTM(m.center, m.radius, [a0], extra + (exact ? 0n : 2n), m.order, P);
}

/**
 * Horner's rule in fixed point at an exact point w with |w| ≤ 1. Returns the
 * value and a bound on its distance from the exact value of the polynomial.
 */
function evalFixed(a, w, P) {
  const wf = fixedOf(w, P);
  let v = a[a.length - 1];
  let S = 0n;
  for (const c of a) S += up(c);
  for (let k = a.length - 2; k >= 0; k--) {
    v = [((v[0] * wf[0] - v[1] * wf[1]) >> P) + a[k][0], ((v[0] * wf[1] + v[1] * wf[0]) >> P) + a[k][1]];
  }
  // each step: the running value (at most S) times the error of wf (under 2 units), plus 2 units of rounding,
  // carried through at most n further multiplications by |wf| ≤ 1 + 2 units
  const n = BigInt(a.length);
  return { v, err: 2n * n * (2n * ((S >> P) + 1n) + 4n) };
}

// ---- functions of a model ----

const pointOf = (b) => b.value ?? b.convention?.value ?? b.candidates[0].value;

function expModel(f) {
  const n = f.order;
  const h = f._tail();
  const B = h._spread();
  const bound = new Rational(BigInt(n + 2));
  if (B.compareTo(bound) >= 0) refuse("precision", `exp: the function varies by up to ${B.toNumber().toPrecision(3)} on ${short(f)}, too much for order ${n}; use a smaller disc or a higher order`);
  let s = f._lift(1);
  for (let k = n; k >= 1; k--) s = h.mul(s).mul(new Rational(1n, BigInt(k))).add(1);
  // Σ_{k>n} B^k/k! ≤ B^(n+1)/(n+1)! · 1/(1 − B/(n+2))
  let fact = 1n;
  for (let k = 2n; k <= BigInt(n + 1); k++) fact *= k;
  const tail = B.pow(n + 1).divide(new Rational(fact)).divide(ONE.subtract(B.divide(bound)));
  return s._more(tail).mul(expPoint(f._c(f.a[0]), f.bits));
}
/** u = (f − f₀)/f₀ and a bound q < 1 for it, for the series of log and of roots around the constant term. */
function relative(f, name) {
  const a0 = f.a[0];
  if (down(a0) === 0n) refuse("branch-point", `${name}: the function may vanish on ${short(f)}`);
  const u = f._tail().mul(f._c(a0).inv());
  const q = new Rational(u._sum() + u.e, 1n << f.P);
  if (q.compareTo(new Rational(7n, 8n)) > 0) {
    refuse(q.compareTo(ONE) >= 0 ? "branch-point" : "precision", `${name}: the function may come too close to 0 on ${short(f)} (it varies by ${q.toNumber().toPrecision(3)} of its central value); use a smaller disc`);
  }
  return { u, q };
}
/**
 * The principal branch is not analytic across the negative real axis, so a
 * model of it on a disc whose image reaches the cut would be a model of
 * something else (the branch continued across), and along a contour the
 * pieces would not even agree with each other. Refused.
 */
function offTheCut(f, name) {
  const [x, y] = f.a[0];
  const B = f._sum(1) + f.e + 2n;
  if (x <= 2n && (y < 0n ? -y : y) <= B) {
    refuse("branch-cut", `${name}: the function may reach the negative real axis on ${short(f)}, where the principal branch jumps. A branch cut elsewhere avoids this: log(-z) + i*pi is a logarithm cut along the positive real axis`);
  }
}
function hornerSeries(f, u, coef) {
  const n = f.order;
  let s = f._lift(coef(n));
  for (let k = n - 1; k >= 0; k--) s = u.mul(s).add(coef(k));
  return s;
}
function logModel(f, opts = {}) {
  const { branch = 0 } = opts;
  const { u, q } = relative(f, "log");
  offTheCut(f, "log");
  const n = f.order;
  // log(1 + u) = Σ_{k≥1} (−1)^(k+1) u^k/k; the tail is at most q^(n+1)/((n+1)(1 − q))
  const s = hornerSeries(f, u, (k) => (k === 0 ? ZERO : new Rational(k % 2 ? 1n : -1n, BigInt(k))));
  const tail = q.pow(n + 1).divide(new Rational(BigInt(n + 1)).multiply(ONE.subtract(q)));
  let L0 = pointOf(logPoint(f._c(f.a[0]), f.bits));
  if (branch) {
    const pi = constants.pi(f.bits);
    const turn = new Rational(BigInt(2 * branch));
    L0 = L0.add(new Rect(iv(ZERO, ZERO), branch > 0 ? iv(pi.low.multiply(turn), pi.high.multiply(turn)) : iv(pi.high.multiply(turn), pi.low.multiply(turn))));
  }
  return s._more(tail).add(L0);
}
function rootModel(f, m, opts = {}) {
  const { branch = 0 } = opts;
  if (!Number.isInteger(m) || m < 1) throw new RangeError("root(f, n) needs a whole number n ≥ 1");
  if (m === 1) return f;
  const { u, q } = relative(f, "root");
  offTheCut(f, "root");
  const n = f.order;
  // (1 + u)^(1/m) = Σ C(1/m, k) u^k with |C(1/m, k)| ≤ 1; the tail is at most q^(n+1)/(1 − q)
  const b = [ONE];
  const alpha = new Rational(1n, BigInt(m));
  for (let k = 1; k <= n; k++) b.push(b[k - 1].multiply(alpha.subtract(new Rational(BigInt(k - 1)))).divide(new Rational(BigInt(k))));
  const s = hornerSeries(f, u, (k) => b[k]);
  const tail = q.pow(n + 1).divide(ONE.subtract(q));
  const r0 = rootPoint(f._c(f.a[0]), m, { bits: f.bits, ...(branch ? { branch } : {}) });
  return s._more(tail).mul(pointOf(r0));
}
function powModel(f, w, opts) {
  if (f instanceof CTM && !(w instanceof CTM)) {
    if (typeof w === "number" && Number.isInteger(w)) return f.pow(w);
    const e = w instanceof Rect || w instanceof Disc ? w : Complex.from(w);
    if (e instanceof Complex && e.isReal() && e.re.denominator === 1n) return f.pow(Number(e.re.numerator));
    return expModel(logModel(f, opts).mul(e));
  }
  const base = f instanceof CTM ? f : w._lift(f);
  return expModel(logModel(base, opts).mul(w));
}
const timesI = (f) => f.mul(Complex.I);
const HALF = new Rational(1n, 2n);
const model = {
  exp: expModel,
  log: logModel,
  sqrt: (f, o) => rootModel(f, 2, o),
  root: rootModel,
  pow: powModel,
  // exp(if) and exp(−if) are built separately: inverting one of them by its geometric series converges far more slowly
  sin(f) {
    const g = timesI(f);
    return expModel(g).sub(expModel(g.neg())).mul(new Complex(ZERO, HALF.negate()));
  },
  cos(f) {
    const g = timesI(f);
    return expModel(g).add(expModel(g.neg())).scale(HALF);
  },
  sinh(f) {
    return expModel(f).sub(expModel(f.neg())).scale(HALF);
  },
  cosh(f) {
    return expModel(f).add(expModel(f.neg())).scale(HALF);
  },
  tan(f) {
    return model.sin(f).div(model.cos(f));
  },
  tanh(f) {
    return model.sinh(f).div(model.cosh(f));
  },
};
/** For special.js: gamma, zeta and their entire companions as functions of models. */
export const _registerModel = (name, fn) => (model[name] = fn);
_setModelHook({
  is: (x) => x instanceof CTM,
  apply(name, args) {
    const fn = model[name];
    if (!fn) refuse("unsupported", `${name} is not available for Taylor models; build it from exp and log`);
    return fn(...args);
  },
});

/**
 * The Taylor model of f on a disc. f is a function of one argument built from
 * the methods add, sub, mul, div, neg, pow, scale and the package's exp, log,
 * sqrt, root, pow, sin, cos, tan, sinh, cosh, tanh; or the same as text,
 * "exp(z) - 3z^2". Options: order (default 12), bits (default 96).
 */
export function taylorModel(f, disc, opts) {
  const z = CTM.z(disc, opts);
  const out = asFunction(f)(z);
  return out instanceof CTM ? out : z._lift(out);
}

/** A function of models from a JavaScript function or from text in the variable z (i, pi and e are understood). */
export function asFunction(f) {
  if (typeof f === "function") return f;
  if (typeof f === "string") return compile(f);
  if (f && typeof f.type === "string") return compileAst(f);
  throw new TypeError("give a function of z, or its formula as text");
}
let parser = null;
/** Set by index.js: analytic's parse, kept out of this file's own imports so the models can be used without it. */
export const _setParser = (p) => (parser = p);
function compile(text) {
  return compileAst(parseFormula(text));
}
/**
 * Text to a syntax tree. The parser knows gamma but not zeta, and would read
 * zeta(z) as a product; zeta goes through it under a borrowed name.
 */
export function parseFormula(text) {
  if (!parser) refuse("unsupported", "formulas as text need analytic's parser");
  if (/\btrigamma\b/.test(text)) refuse("unsupported", "trigamma is not available here");
  const swap = (node) => (node.type === "call" && node.name === "trigamma" ? { ...node, name: "zeta", args: node.args.map(swap) } : node.args ? { ...node, args: node.args.map(swap) } : node);
  return swap(parser(text.replace(/\bzeta\s*\(/g, "trigamma(")));
}
function compileAst(ast) {
  const names = { ln: "log" };
  const run = (node, z) => {
    switch (node.type) {
      case "num":
        return z._lift(node.value);
      case "var":
        if (node.name === "z") return z;
        if (node.name === "i") return z._lift(Complex.I);
        return refuse("bad-argument", `unknown name ${node.name}: the variable is z`);
      case "const": {
        const c = node.name === "pi" ? constants.pi(z.bits) : constants.e(z.bits);
        return z._lift(new Rect(iv(c.low, c.high), iv(ZERO, ZERO)));
      }
      case "op": {
        const args = node.args.map((a) => run(a, z));
        switch (node.op) {
          case "+": return args[0].add(args[1]);
          case "-": return args[0].sub(args[1]);
          case "*": return args[0].mul(args[1]);
          case "/": return args[0].div(args[1]);
          case "neg": return args[0].neg();
          case "^": {
            const e = node.args[1];
            if (e.type === "num" && e.value.denominator === 1n) return args[0].pow(Number(e.value.numerator));
            if (e.type === "op" && e.op === "neg" && e.args[0].type === "num" && e.args[0].value.denominator === 1n) return args[0].pow(-Number(e.args[0].value.numerator));
            return powModel(args[0], args[1]);
          }
        }
        return refuse("unsupported", `operator ${node.op}`);
      }
      case "call": {
        const name = names[node.name] ?? node.name;
        if (!model[name] || name === "root" || name === "pow") return refuse("unsupported", `${node.name} is not available here; the functions are exp, log, sqrt, sin, cos, tan, sinh, cosh, tanh, gamma, zeta`);
        return model[name](run(node.args[0], z));
      }
    }
    return refuse("unsupported", `cannot read a ${node.type} node`);
  };
  return (z) => run(ast, z);
}
