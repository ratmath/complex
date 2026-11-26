// Exact complex rationals: a + bi with a and b Core Rationals.
//
// Every operation here is exact. Numbers with square roots or transcendental
// parts (a modulus, an argument, exp of anything) are enclosures instead:
// see enclosure.js and functions.js.

import { Rational } from "@ratmath/core";

export class ComplexError extends Error {
  constructor(code, message) {
    super(message);
    this.name = "ComplexError";
    this.code = code;
  }
}
export const refuse = (code, message) => {
  throw new ComplexError(code, message);
};

export const ZERO = new Rational(0n);
export const ONE = new Rational(1n);
export const TWO = new Rational(2n);

/** A Core Rational from a Rational, bigint, whole or decimal number, or text in Core's notation ("1/2", "0.#3", "1..2/3"). */
export function Q(x) {
  if (x instanceof Rational) return x;
  if (typeof x === "bigint") return new Rational(x);
  if (typeof x === "number") {
    if (!Number.isFinite(x)) throw new RangeError(`not a finite number: ${x}`);
    // the shortest decimal that names the double: 0.1 is 1/10
    const [mant, exp = "0"] = String(x).split("e");
    const q = new Rational(mant);
    const e = Number(exp);
    return e === 0 ? q : e > 0 ? q.multiply(new Rational(10n ** BigInt(e))) : q.divide(new Rational(10n ** BigInt(-e)));
  }
  if (typeof x === "string") return new Rational(x.trim());
  if (x && typeof x.numerator === "bigint") return new Rational(x.numerator, x.denominator);
  throw new TypeError(`cannot read ${x} as a rational number`);
}

const fmt = (q) => q.toString();
const sq = (q) => q.multiply(q);

/** floor(sqrt(n)) for a non-negative bigint. */
export function isqrt(n) {
  if (n < 2n) return n;
  let x = 1n << BigInt((n.toString(2).length + 1) >> 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
}

/** The exact square root of a non-negative rational when it is rational, else null. */
export function exactSqrt(q) {
  if (q.numerator < 0n) return null;
  const a = isqrt(q.numerator);
  const b = isqrt(q.denominator);
  return a * a !== q.numerator || b * b !== q.denominator ? null : new Rational(a, b);
}

export class Complex {
  constructor(re = ZERO, im = ZERO) {
    this.re = Q(re);
    this.im = Q(im);
    Object.freeze(this);
  }

  /** From a Complex, a real number in any form Q reads, a [re, im] pair, or text such as "3 - 4i" or "(1 + 2i)/(3 - i)". */
  static from(x) {
    if (x instanceof Complex) return x;
    if (Array.isArray(x) && x.length === 2) return new Complex(x[0], x[1]);
    if (typeof x === "string") return Complex.parse(x);
    if (x && typeof x === "object" && "re" in x && "im" in x) return new Complex(x.re, x.im);
    return new Complex(Q(x), ZERO);
  }

  /**
   * Read exact complex arithmetic: numbers in Core's notation, `i`, + - * /,
   * ^ with a whole-number exponent, and parentheses. "2i", "1/2 - 3/4i",
   * "(1 + i)^8", "0.#3 + 1..1/2 i".
   */
  static parse(text) {
    const toks = [];
    const re = /\s*(?:(\d+\.\.\d+\/\d+|\d*\.?\d*#\d+|\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+)|([ij](?![A-Za-z]))|(\*\*|[-+*/^()]))/y;
    let pos = 0;
    while (pos < text.length) {
      if (/^\s*$/.test(text.slice(pos))) break;
      re.lastIndex = pos;
      const m = re.exec(text);
      if (!m) throw new SyntaxError(`cannot read ${JSON.stringify(text.slice(pos))} in ${JSON.stringify(text)}`);
      toks.push(m[1] !== undefined ? { t: "num", v: m[1] } : m[2] ? { t: "i" } : { t: m[3] === "**" ? "^" : m[3] });
      pos = re.lastIndex;
    }
    let k = 0;
    const peek = () => toks[k]?.t;
    const expr = () => {
      let left = term();
      while (peek() === "+" || peek() === "-") left = toks[k++].t === "+" ? left.add(term()) : left.sub(term());
      return left;
    };
    const starts = () => peek() === "num" || peek() === "i" || peek() === "(";
    const term = () => {
      let left = unary();
      for (;;) {
        if (peek() === "*") {
          k++;
          left = left.mul(unary());
        } else if (peek() === "/") {
          k++;
          // "3/4i" reads as (3/4)·i, the way it is written by hand
          const d = power();
          left = left.div(d);
          if (peek() === "i" ) {
            k++;
            left = left.mul(Complex.I);
          }
        } else if (starts()) left = left.mul(power());
        else return left;
      }
    };
    const unary = () => {
      if (peek() === "-") {
        k++;
        return unary().neg();
      }
      if (peek() === "+") {
        k++;
        return unary();
      }
      return power();
    };
    const power = () => {
      const base = primary();
      if (peek() !== "^") return base;
      k++;
      let sign = 1;
      if (peek() === "-") {
        k++;
        sign = -1;
      }
      const e = toks[k++];
      if (!e || e.t !== "num" || !/^\d+$/.test(e.v)) throw new SyntaxError(`a power must be a whole number in ${JSON.stringify(text)}`);
      return base.pow(sign * Number(e.v));
    };
    const primary = () => {
      const tok = toks[k++];
      if (!tok) throw new SyntaxError(`unexpected end of ${JSON.stringify(text)}`);
      if (tok.t === "num") return new Complex(new Rational(tok.v), ZERO);
      if (tok.t === "i") return Complex.I;
      if (tok.t === "(") {
        const e = expr();
        if (toks[k++]?.t !== ")") throw new SyntaxError(`missing ) in ${JSON.stringify(text)}`);
        return e;
      }
      throw new SyntaxError(`unexpected ${tok.t} in ${JSON.stringify(text)}`);
    };
    const out = expr();
    if (k < toks.length) throw new SyntaxError(`unexpected ${toks[k].t} in ${JSON.stringify(text)}`);
    return out;
  }

  add(o) {
    o = Complex.from(o);
    return new Complex(this.re.add(o.re), this.im.add(o.im));
  }
  sub(o) {
    o = Complex.from(o);
    return new Complex(this.re.subtract(o.re), this.im.subtract(o.im));
  }
  mul(o) {
    o = Complex.from(o);
    return new Complex(this.re.multiply(o.re).subtract(this.im.multiply(o.im)), this.re.multiply(o.im).add(this.im.multiply(o.re)));
  }
  div(o) {
    o = Complex.from(o);
    const n = o.norm();
    if (n.equals(ZERO)) refuse("division-by-zero", "division by zero");
    const p = this.mul(o.conj());
    return new Complex(p.re.divide(n), p.im.divide(n));
  }
  neg() {
    return new Complex(this.re.negate(), this.im.negate());
  }
  /** The conjugate a − bi. */
  conj() {
    return new Complex(this.re, this.im.negate());
  }
  inv() {
    return Complex.ONE.div(this);
  }
  /** |z|² = a² + b², exact. (The modulus itself is usually irrational: see modulus().) */
  norm() {
    return sq(this.re).add(sq(this.im));
  }
  /** |z| when it is rational (3 + 4i has modulus 5), else null. */
  exactModulus() {
    return exactSqrt(this.norm());
  }
  /** A whole-number power, exact. */
  pow(n) {
    if (typeof n === "bigint") n = Number(n);
    if (!Number.isInteger(n)) throw new RangeError("Complex.pow takes a whole-number exponent; use pow() from the functions for others");
    if (n < 0) return this.inv().pow(-n);
    let result = Complex.ONE;
    let base = this;
    while (n > 0) {
      if (n & 1) result = result.mul(base);
      n >>= 1;
      if (n) base = base.mul(base);
    }
    return result;
  }
  /** Multiply by a real scalar. */
  scale(q) {
    q = Q(q);
    return new Complex(this.re.multiply(q), this.im.multiply(q));
  }
  equals(o) {
    o = Complex.from(o);
    return this.re.equals(o.re) && this.im.equals(o.im);
  }
  isZero() {
    return this.re.equals(ZERO) && this.im.equals(ZERO);
  }
  isReal() {
    return this.im.equals(ZERO);
  }
  isImaginary() {
    return this.re.equals(ZERO) && !this.im.equals(ZERO);
  }
  /** Both parts whole numbers. */
  isGaussianInteger() {
    return this.re.denominator === 1n && this.im.denominator === 1n;
  }
  /**
   * The exact square roots when they are complex rationals: [w, −w] with
   * Re w ≥ 0 (the principal one first), else null. sqrt(3 + 4i) = 2 + i.
   */
  exactSqrt() {
    if (this.isZero()) return [Complex.ZERO, Complex.ZERO];
    const m = this.exactModulus();
    if (m === null) return null;
    // w = x + iy with x² = (m + a)/2, y² = (m − a)/2, sign(y) = sign(b)
    const x = exactSqrt(m.add(this.re).divide(TWO));
    const y = exactSqrt(m.subtract(this.re).divide(TWO));
    if (x === null || y === null) return null;
    const w = new Complex(x, this.im.numerator < 0n ? y.negate() : y);
    return [w, w.neg()];
  }
  /** [re, im] as doubles. */
  toNumbers() {
    return [this.re.toNumber(), this.im.toNumber()];
  }
  toString() {
    const { re, im } = this;
    if (im.equals(ZERO)) return fmt(re);
    const a = im.abs();
    const part = a.equals(ONE) ? "i" : `${fmt(a)}i`;
    if (re.equals(ZERO)) return im.numerator < 0n ? `-${part}` : part;
    return `${fmt(re)} ${im.numerator < 0n ? "-" : "+"} ${part}`;
  }
  toTeX() {
    const t = (q) => (q.denominator === 1n ? `${q.numerator}` : `${q.numerator < 0n ? "-" : ""}\\frac{${q.numerator < 0n ? -q.numerator : q.numerator}}{${q.denominator}}`);
    const { re, im } = this;
    if (im.equals(ZERO)) return t(re);
    const a = im.abs();
    const part = a.equals(ONE) ? "i" : `${t(a)}i`;
    if (re.equals(ZERO)) return im.numerator < 0n ? `-${part}` : part;
    return `${t(re)} ${im.numerator < 0n ? "-" : "+"} ${part}`;
  }
  toJSON() {
    return { re: this.re.toString(), im: this.im.toString() };
  }
}
Complex.ZERO = new Complex(ZERO, ZERO);
Complex.ONE = new Complex(ONE, ZERO);
Complex.I = new Complex(ZERO, ONE);

/** A complex rational: complex("3 - 4i"), complex(3, -4), complex([3, -4]). */
export const complex = (re, im) => (im === undefined ? Complex.from(re) : new Complex(re, im));

/**
 * The complex rationals Q(i) as a RatMath field object, for linalg matrices
 * and structures polynomials with complex entries.
 */
export const CQ = Object.freeze({
  name: "Q(i)",
  zero: Complex.ZERO,
  one: Complex.ONE,
  add: (a, b) => Complex.from(a).add(b),
  sub: (a, b) => Complex.from(a).sub(b),
  mul: (a, b) => Complex.from(a).mul(b),
  div: (a, b) => Complex.from(a).div(b),
  neg: (a) => Complex.from(a).neg(),
  inv: (a) => Complex.from(a).inv(),
  /** 0 for zero, 1 otherwise: the complex numbers have no order. */
  sign: (a) => (Complex.from(a).isZero() ? 0 : 1),
  equals: (a, b) => Complex.from(a).equals(b),
  from: (x) => Complex.from(x),
  of: (x) => Complex.from(x),
  conj: (a) => Complex.from(a).conj(),
  format: (a) => Complex.from(a).toString(),
  characteristic: 0,
});
