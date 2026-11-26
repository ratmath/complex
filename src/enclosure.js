// Enclosures of complex numbers with rational data: rectangles and discs.
//
// An enclosure is a set that is certain to contain the true value. Arithmetic
// on enclosures returns an enclosure of every possible result, so whatever is
// computed from them stays certain.

import { Rational, RationalInterval } from "@ratmath/core";
import { Complex, Q, ZERO, ONE, TWO, isqrt, refuse } from "./complex.js";


// ---- real intervals as { lo, hi } --------------------------------------------------

export const iv = (lo, hi = lo) => ({ lo, hi });
const min = (a, b) => (a.compareTo(b) <= 0 ? a : b);
const max = (a, b) => (a.compareTo(b) >= 0 ? a : b);
export const I = {
  of(x) {
    if (x instanceof RationalInterval) return iv(x.low, x.high);
    if (x && x.lo instanceof Rational) return x;
    if (Array.isArray(x)) return iv(min(Q(x[0]), Q(x[1])), max(Q(x[0]), Q(x[1])));
    const q = Q(x);
    return iv(q, q);
  },
  add: (a, b) => iv(a.lo.add(b.lo), a.hi.add(b.hi)),
  sub: (a, b) => iv(a.lo.subtract(b.hi), a.hi.subtract(b.lo)),
  neg: (a) => iv(a.hi.negate(), a.lo.negate()),
  mul(a, b) {
    const p = [a.lo.multiply(b.lo), a.lo.multiply(b.hi), a.hi.multiply(b.lo), a.hi.multiply(b.hi)];
    return iv(p.reduce(min), p.reduce(max));
  },
  scale: (a, q) => (q.numerator >= 0n ? iv(a.lo.multiply(q), a.hi.multiply(q)) : iv(a.hi.multiply(q), a.lo.multiply(q))),
  sq(a) {
    const l = a.lo.multiply(a.lo);
    const h = a.hi.multiply(a.hi);
    return I.hasZero(a) ? iv(ZERO, max(l, h)) : iv(min(l, h), max(l, h));
  },
  recip(a) {
    if (I.hasZero(a)) refuse("division-by-zero", "division by an interval that contains zero");
    return iv(ONE.divide(a.hi), ONE.divide(a.lo));
  },
  div: (a, b) => I.mul(a, I.recip(b)),
  hull: (a, b) => iv(min(a.lo, b.lo), max(a.hi, b.hi)),
  meet(a, b) {
    const lo = max(a.lo, b.lo);
    const hi = min(a.hi, b.hi);
    return lo.compareTo(hi) <= 0 ? iv(lo, hi) : null;
  },
  hasZero: (a) => a.lo.numerator <= 0n && a.hi.numerator >= 0n,
  has: (a, q) => a.lo.compareTo(q) <= 0 && a.hi.compareTo(q) >= 0,
  width: (a) => a.hi.subtract(a.lo),
  isPoint: (a) => a.lo.equals(a.hi),
  mid: (a) => a.lo.add(a.hi).divide(TWO),
  /** The same interval with its ends moved outward to fractions with about `bits` significant bits. */
  round: (a, bits) => iv(roundDown(a.lo, bits), roundUp(a.hi, bits)),
  toCore: (a) => new RationalInterval(a.lo, a.hi),
};

const bitlen = (n) => (n < 0n ? -n : n).toString(2).length;
/** The largest fraction ≤ q with a power-of-two denominator and about `bits` significant bits. */
export function roundDown(q, bits) {
  if (q.numerator === 0n || bitlen(q.denominator) <= bits) return q;
  const shift = BigInt(Math.max(0, bits + bitlen(q.denominator) - bitlen(q.numerator) + 1));
  const scaled = q.numerator << shift;
  let f = scaled / q.denominator;
  if (scaled % q.denominator !== 0n && scaled < 0n) f -= 1n;
  return new Rational(f, 1n << shift);
}
export const roundUp = (q, bits) => roundDown(q.negate(), bits).negate();

/** An upper bound for sqrt(q), q ≥ 0, good to about `bits` bits; exact when the root is rational. */
export function sqrtUp(q, bits = 64) {
  if (q.numerator === 0n) return ZERO;
  const shift = BigInt(2 * Math.max(8, bits + bitlen(q.denominator) - bitlen(q.numerator)));
  const n = (q.numerator << shift) / q.denominator; // floor(q·2^shift)
  const r = isqrt(n);
  const exact = r * r === n && (q.numerator << shift) % q.denominator === 0n;
  return new Rational(exact ? r : r + 1n, 1n << (shift / 2n));
}
/** A lower bound for sqrt(q). */
export function sqrtDown(q, bits = 64) {
  if (q.numerator === 0n) return ZERO;
  const shift = BigInt(2 * Math.max(8, bits + bitlen(q.denominator) - bitlen(q.numerator)));
  return new Rational(isqrt((q.numerator << shift) / q.denominator), 1n << (shift / 2n));
}

/** q as a decimal rounded to `places` places after the point, without trailing zeros. */
const fixed = (q, places) => {
  const neg = q.numerator < 0n;
  const a = q.abs();
  const scaled = (a.numerator * 10n ** BigInt(places) * 2n + a.denominator) / (a.denominator * 2n);
  let t = scaled.toString().padStart(places + 1, "0");
  t = places ? `${t.slice(0, t.length - places)}.${t.slice(t.length - places)}`.replace(/\.?0+$/, "") : t;
  return (neg && scaled !== 0n ? "-" : "") + (t || "0");
};
/** How many places after the point give `digits` significant digits of the larger of the given numbers. */
const placesFor = (digits, ...qs) => {
  let e = -Infinity;
  for (const q of qs) if (q.numerator !== 0n) e = Math.max(e, q.abs().numerator.toString().length - q.abs().denominator.toString().length);
  return Number.isFinite(e) ? Math.max(0, digits - e - 1) : 0;
};
const dec = (q, digits) => (q.denominator === 1n ? q.numerator.toString() : fixed(q, placesFor(digits, q)));

// ---- rectangles -----------------------------------------------------------------------

/** The rectangle [re.lo, re.hi] + i[im.lo, im.hi] with rational sides. */
export class Rect {
  constructor(re, im) {
    this.re = I.of(re);
    this.im = I.of(im);
    Object.freeze(this.re);
    Object.freeze(this.im);
    Object.freeze(this);
  }
  /** From a Rect, a Disc (its bounding square), a Complex or anything Complex.from reads. */
  static from(x) {
    if (x instanceof Rect) return x;
    if (x instanceof Disc) return x.toRect();
    const z = Complex.from(x);
    return new Rect(iv(z.re), iv(z.im));
  }
  /** The smallest rectangle around the given points or enclosures. */
  static hull(...items) {
    return items.map(Rect.from).reduce((a, b) => a.hull(b));
  }
  get isPoint() {
    return I.isPoint(this.re) && I.isPoint(this.im);
  }
  /** The midpoint, an exact Complex. */
  get center() {
    return new Complex(I.mid(this.re), I.mid(this.im));
  }
  get width() {
    return I.width(this.re);
  }
  get height() {
    return I.width(this.im);
  }
  /** The larger of width and height. */
  get size() {
    return this.width.compareTo(this.height) >= 0 ? this.width : this.height;
  }
  /** The four corners, counter-clockwise from the lower left. */
  corners() {
    const { re, im } = this;
    return [new Complex(re.lo, im.lo), new Complex(re.hi, im.lo), new Complex(re.hi, im.hi), new Complex(re.lo, im.hi)];
  }
  add(o) {
    o = Rect.from(o);
    return new Rect(I.add(this.re, o.re), I.add(this.im, o.im));
  }
  sub(o) {
    o = Rect.from(o);
    return new Rect(I.sub(this.re, o.re), I.sub(this.im, o.im));
  }
  mul(o) {
    o = Rect.from(o);
    return new Rect(I.sub(I.mul(this.re, o.re), I.mul(this.im, o.im)), I.add(I.mul(this.re, o.im), I.mul(this.im, o.re)));
  }
  neg() {
    return new Rect(I.neg(this.re), I.neg(this.im));
  }
  conj() {
    return new Rect(this.re, I.neg(this.im));
  }
  /** Multiply by a real scalar. */
  scale(q) {
    q = Q(q);
    return new Rect(I.scale(this.re, q), I.scale(this.im, q));
  }
  /** Multiply by i: a quarter turn. */
  timesI() {
    return new Rect(I.neg(this.im), this.re);
  }
  /** The range of |z|² over the rectangle, as { lo, hi }. */
  normRange() {
    return I.add(I.sq(this.re), I.sq(this.im));
  }
  inv() {
    if (this.containsZero()) refuse("division-by-zero", `division by an enclosure that contains zero: ${this}`);
    const n = this.normRange();
    return new Rect(I.div(this.re, n), I.div(I.neg(this.im), n));
  }
  div(o) {
    o = Rect.from(o);
    if (o.isPoint && this.isPoint) return Rect.from(this.center.div(o.center));
    return this.mul(o.inv());
  }
  /** A whole-number power by repeated multiplication. */
  pow(n) {
    if (!Number.isInteger(n)) throw new RangeError("Rect.pow takes a whole-number exponent");
    if (n < 0) return this.inv().pow(-n);
    let result = Rect.from(Complex.ONE);
    let base = this;
    while (n > 0) {
      if (n & 1) result = result.mul(base);
      n >>= 1;
      if (n) base = base.mul(base);
    }
    return result;
  }
  containsZero() {
    return I.hasZero(this.re) && I.hasZero(this.im);
  }
  /** True when the point or the whole enclosure lies inside. */
  contains(o) {
    o = Rect.from(o);
    return this.re.lo.compareTo(o.re.lo) <= 0 && this.re.hi.compareTo(o.re.hi) >= 0 && this.im.lo.compareTo(o.im.lo) <= 0 && this.im.hi.compareTo(o.im.hi) >= 0;
  }
  intersects(o) {
    return this.intersect(o) !== null;
  }
  /** The common part, or null. */
  intersect(o) {
    o = Rect.from(o);
    const re = I.meet(this.re, o.re);
    const im = I.meet(this.im, o.im);
    return re && im ? new Rect(re, im) : null;
  }
  hull(o) {
    o = Rect.from(o);
    return new Rect(I.hull(this.re, o.re), I.hull(this.im, o.im));
  }
  /** The same set, or slightly more, with shorter fractions (about `bits` significant bits). */
  round(bits = 64) {
    return new Rect(I.round(this.re, bits), I.round(this.im, bits));
  }
  /** The disc through the corners, centred at the midpoint. */
  toDisc() {
    const half = new Complex(this.width.divide(TWO), this.height.divide(TWO));
    return new Disc(this.center, sqrtUp(half.norm()));
  }
  /** Four quarters. */
  split() {
    const c = this.center;
    const { re, im } = this;
    return [new Rect(iv(re.lo, c.re), iv(im.lo, c.im)), new Rect(iv(c.re, re.hi), iv(im.lo, c.im)), new Rect(iv(c.re, re.hi), iv(c.im, im.hi)), new Rect(iv(re.lo, c.re), iv(c.im, im.hi))];
  }
  /** Real and imaginary parts as Core RationalIntervals. */
  get real() {
    return I.toCore(this.re);
  }
  get imag() {
    return I.toCore(this.im);
  }
  /** Decimal text. A point prints exactly; otherwise the digits shown are those of the midpoint, with the half-width after ±. */
  toString(digits = 12) {
    if (this.isPoint) return this.center.toString();
    const c = this.center;
    const r = this.size.divide(TWO);
    const tiny = c.re.abs().compareTo(c.im.abs()) > 0 ? c.re.abs() : c.im.abs();
    // a rectangle about as close to 0 as it is wide: written around 0, with the midpoint's offset added to the ±
    if (tiny.compareTo(r) <= 0) return `0 + 0i ± ${Number(r.add(tiny).toNumber().toPrecision(2)).toExponential()}`;
    // both parts to the same number of places, set by the larger one: a part that is 0 to that accuracy prints as 0
    const places = placesFor(digits, c.re, c.im);
    const re = fixed(c.re, places);
    const im = fixed(c.im.abs(), places);
    const sign = c.im.numerator < 0n && im !== "0" ? "-" : "+";
    return `${re} ${sign} ${im}i ± ${Number(r.toNumber().toPrecision(2)).toExponential()}`;
  }
  toJSON() {
    return { re: [this.re.lo.toString(), this.re.hi.toString()], im: [this.im.lo.toString(), this.im.hi.toString()] };
  }
}

/** A rectangle: rect("1 + 2i"), rect([0, 1], [2, 3]) for [0,1] + i[2,3], rect(z, radius) for the square around z. */
export function rect(a, b) {
  if (b === undefined) return Rect.from(a);
  if (Array.isArray(a) || a instanceof RationalInterval) return new Rect(I.of(a), I.of(b));
  const z = Complex.from(a);
  const r = Q(b);
  return new Rect(iv(z.re.subtract(r), z.re.add(r)), iv(z.im.subtract(r), z.im.add(r)));
}

// ---- discs -------------------------------------------------------------------------

/** The closed disc |z − center| ≤ radius with rational centre and radius. */
export class Disc {
  constructor(center, radius = ZERO) {
    this.center = Complex.from(center);
    this.radius = Q(radius);
    if (this.radius.numerator < 0n) throw new RangeError("a disc needs a radius ≥ 0");
    Object.freeze(this);
  }
  static from(x) {
    if (x instanceof Disc) return x;
    if (x instanceof Rect) return x.toDisc();
    return new Disc(Complex.from(x), ZERO);
  }
  get isPoint() {
    return this.radius.equals(ZERO);
  }
  add(o) {
    o = Disc.from(o);
    return new Disc(this.center.add(o.center), this.radius.add(o.radius));
  }
  sub(o) {
    o = Disc.from(o);
    return new Disc(this.center.sub(o.center), this.radius.add(o.radius));
  }
  neg() {
    return new Disc(this.center.neg(), this.radius);
  }
  conj() {
    return new Disc(this.center.conj(), this.radius);
  }
  /** |c₁c₂ − z₁z₂| ≤ |c₁|r₂ + |c₂|r₁ + r₁r₂. */
  mul(o) {
    o = Disc.from(o);
    const m1 = this.isPoint && o.isPoint ? ZERO : sqrtUp(this.center.norm());
    const m2 = this.isPoint && o.isPoint ? ZERO : sqrtUp(o.center.norm());
    return new Disc(this.center.mul(o.center), m1.multiply(o.radius).add(m2.multiply(this.radius)).add(this.radius.multiply(o.radius)));
  }
  /** The image under z → 1/z, which is again a disc, exactly: centre c̄/(|c|² − r²), radius r/(|c|² − r²). */
  inv() {
    const d = this.center.norm().subtract(this.radius.multiply(this.radius));
    if (d.numerator <= 0n) refuse("division-by-zero", `division by a disc that contains zero: ${this}`);
    return new Disc(this.center.conj().scale(ONE.divide(d)), this.radius.divide(d));
  }
  div(o) {
    return this.mul(Disc.from(o).inv());
  }
  scale(q) {
    q = Q(q);
    return new Disc(this.center.scale(q), this.radius.multiply(q.abs()));
  }
  /** Exact: a point is inside when |z − c|² ≤ r²; a disc when |c₁ − c₂| + r₂ ≤ r₁. */
  contains(o) {
    if (o instanceof Rect) return o.corners().every((z) => this.contains(z));
    o = Disc.from(o);
    const gap = this.radius.subtract(o.radius);
    return gap.numerator >= 0n && this.center.sub(o.center).norm().compareTo(gap.multiply(gap)) <= 0;
  }
  containsZero() {
    return this.contains(Complex.ZERO);
  }
  /** Exact: two discs meet when |c₁ − c₂| ≤ r₁ + r₂. */
  intersects(o) {
    o = Disc.from(o);
    const sum = this.radius.add(o.radius);
    return this.center.sub(o.center).norm().compareTo(sum.multiply(sum)) <= 0;
  }
  /** A disc containing both. */
  hull(o) {
    o = Disc.from(o);
    if (this.contains(o)) return this;
    if (o.contains(this)) return o;
    // centre c₁ + t(c₂ − c₁); any t in [0, 1] is sound once the radius covers both discs
    const d = sqrtUp(this.center.sub(o.center).norm());
    const t = d.add(o.radius).subtract(this.radius).divide(d.multiply(TWO));
    const r1 = t.multiply(d).add(this.radius);
    const r2 = ONE.subtract(t).multiply(d).add(o.radius);
    return new Disc(this.center.add(o.center.sub(this.center).scale(t)), r1.compareTo(r2) >= 0 ? r1 : r2);
  }
  /** The bounding square. */
  toRect() {
    const { center: c, radius: r } = this;
    return new Rect(iv(c.re.subtract(r), c.re.add(r)), iv(c.im.subtract(r), c.im.add(r)));
  }
  toString(digits = 12) {
    if (this.isPoint) return this.center.toString();
    return `disc(${this.center}, ${dec(this.radius, digits)})`;
  }
  toJSON() {
    return { center: this.center.toJSON(), radius: this.radius.toString() };
  }
}

/** A disc: disc("1 + i", "1/10"). */
export const disc = (center, radius = ZERO) => new Disc(center, radius);
