// Cayley–Dickson numbers: the doubling construction that builds the complex
// numbers from the reals, the quaternions from the complex numbers, the
// octonions from the quaternions, the sedenions from the octonions, and so on.
//
// A pair (a, b) of numbers of one level is a number of the next, with
//
//   (a, b)·(c, d) = (a·c − d*·b,  d·a + b·c*)      (a, b)* = (a*, −b)
//
// Each doubling loses a property: the complex numbers are not ordered, the
// quaternions do not commute, the octonions do not associate, and the
// sedenions have zero divisors. Everything here is exact, over the rationals
// by default or over any RatMath field object.

import { Rational } from "@ratmath/core";
import { Q, ZERO, ONE, refuse } from "./complex.js";

/** The rationals as a field object over Core Rationals. */
export const QQ = Object.freeze({
  name: "Q",
  zero: ZERO,
  one: ONE,
  add: (a, b) => a.add(b),
  sub: (a, b) => a.subtract(b),
  mul: (a, b) => a.multiply(b),
  div: (a, b) => a.divide(b),
  neg: (a) => a.negate(),
  sign: (a) => (a.numerator > 0n ? 1 : a.numerator < 0n ? -1 : 0),
  from: (x) => Q(x),
  format: (a) => a.toString(),
});

/** The multiplication table of the basis: table[i][j] = [index, sign] with e_i·e_j = sign·e_index. */
function basisTable(dim) {
  // doubled one level at a time: with e_i = (e_i, 0) and e_{h+i} = (0, e_i),
  //   (a, 0)(c, 0) = (ac, 0)      (a, 0)(0, d) = (0, da)
  //   (0, b)(c, 0) = (0, b·c*)    (0, b)(0, d) = (−d*·b, 0)
  let T = [[[0, 1]]];
  for (let h = 1; h < dim; h *= 2) {
    const N = Array.from({ length: 2 * h }, () => Array(2 * h));
    const cs = (j) => (j === 0 ? 1 : -1);
    for (let i = 0; i < h; i++) for (let j = 0; j < h; j++) {
      N[i][j] = T[i][j];
      N[i][h + j] = [h + T[j][i][0], T[j][i][1]];
      N[h + i][j] = [h + T[i][j][0], cs(j) * T[i][j][1]];
      N[h + i][h + j] = [T[j][i][0], -cs(j) * T[j][i][1]];
    }
    T = N;
  }
  return T;
}

const NAMES = { 1: "real numbers", 2: "complex numbers", 4: "quaternions", 8: "octonions", 16: "sedenions" };
const LABELS = { 1: ["1"], 2: ["1", "i"], 4: ["1", "i", "j", "k"] };

/** One level of the construction: an algebra of dimension 2^level over a field. */
export class CayleyDickson {
  constructor(level, field = QQ) {
    if (!Number.isInteger(level) || level < 0 || level > 8) throw new RangeError("the level must be a whole number from 0 (reals) to 8 (dimension 256)");
    this.level = level;
    this.dim = 2 ** level;
    this.field = field;
    this.name = NAMES[this.dim] ?? `Cayley–Dickson algebra of dimension ${this.dim}`;
    this.labels = LABELS[this.dim] ?? Array.from({ length: this.dim }, (_, i) => (i === 0 ? "1" : `e${i}`));
    this.table = basisTable(this.dim);
    this.zero = new Hyper(this, Array(this.dim).fill(field.zero));
    this.one = this.basis(0);
    Object.freeze(this);
  }
  /** The i-th basis element: basis(0) is 1. */
  basis(i) {
    const { field: K } = this;
    return new Hyper(this, Array.from({ length: this.dim }, (_, k) => (k === i ? K.one : K.zero)));
  }
  /** All basis elements. */
  get units() {
    return Array.from({ length: this.dim }, (_, i) => this.basis(i));
  }
  /** An element from its coefficients (missing ones are 0), from a scalar, from an element of a lower level, or (over the rationals) from text as parse reads it. */
  of(...coeffs) {
    if (coeffs.length === 1 && coeffs[0] instanceof Hyper) {
      const h = coeffs[0];
      if (h.algebra === this) return h;
      if (h.algebra.dim > this.dim) throw new RangeError(`cannot put an element of the ${h.algebra.name} into the ${this.name}`);
      coeffs = h.coeffs;
    } else if (coeffs.length === 1 && Array.isArray(coeffs[0])) coeffs = coeffs[0];
    else if (coeffs.length === 1 && typeof coeffs[0] === "string" && /[a-z]/i.test(coeffs[0]) && this.field === QQ) return this.parse(coeffs[0]);
    if (coeffs.length > this.dim) throw new RangeError(`the ${this.name} have ${this.dim} coefficients, not ${coeffs.length}`);
    const K = this.field;
    const read = K.from ?? ((x) => x);
    return new Hyper(this, Array.from({ length: this.dim }, (_, k) => (k < coeffs.length ? read(coeffs[k]) : K.zero)));
  }
  /** Read "1 + 2i - 3/4 k" (quaternions) or "e1 + 2 e7" (any level). Rational coefficients only. */
  parse(text) {
    const coeffs = Array(this.dim).fill(null).map(() => ZERO);
    const src = text.replace(/\s+/g, "");
    const re = /([+-]?)(\d+(?:\/\d+)?|\d*\.\d+)?\*?(e\d+|[ijk1])?/gy;
    let pos = 0;
    while (pos < src.length) {
      re.lastIndex = pos;
      const m = re.exec(src);
      if (!m || m[0] === "" || (m[2] === undefined && m[3] === undefined)) throw new SyntaxError(`cannot read ${JSON.stringify(text)} as an element of the ${this.name}`);
      let label = m[3] ?? "1";
      let c = m[2] !== undefined ? Q(m[2]) : ONE;
      // "1" alone is the number, not a label after a missing coefficient
      if (label === "1" && m[2] === undefined) c = ONE;
      const idx = label === "1" ? 0 : label[0] === "e" ? Number(label.slice(1)) : "1ijk".indexOf(label);
      if (idx < 0 || idx >= this.dim) throw new SyntaxError(`${label} is not a unit of the ${this.name}`);
      coeffs[idx] = coeffs[idx].add(m[1] === "-" ? c.negate() : c);
      pos = re.lastIndex;
    }
    return this.of(coeffs);
  }
  /**
   * Which laws hold, found by exact computation on the basis elements and
   * their pairwise sums: { commutative, associative, alternative, flexible,
   * composition, zeroDivisors } where zeroDivisors is null or a pair [a, b]
   * of nonzero elements with a·b = 0. (Commutativity and associativity are
   * linear in each argument, so the basis decides them. The others are checked
   * on basis elements and sums of two, which is where they first fail.)
   */
  properties() {
    const { dim, table: T } = this;
    // sparse integer vectors as Maps index → coefficient
    const mul = (a, b) => {
      const out = new Map();
      for (const [i, x] of a) for (const [j, y] of b) {
        const [k, s] = T[i][j];
        const val = (out.get(k) ?? 0) + s * x * y;
        if (val) out.set(k, val);
        else out.delete(k);
      }
      return out;
    };
    const same = (a, b) => a.size === b.size && [...a].every(([k, x]) => b.get(k) === x);
    const e = (i) => new Map([[i, 1]]);
    const basis = Array.from({ length: dim }, (_, i) => e(i));
    const lim = Math.min(dim, 16);
    const some = basis.slice(0, lim);
    const pairs = [];
    for (let i = 0; i < lim; i++) for (let j = i + 1; j < lim; j++) pairs.push(new Map([[i, 1], [j, 1]]), new Map([[i, 1], [j, -1]]));
    const commutative = basis.every((a, i) => basis.every((b, j) => T[i][j][0] === T[j][i][0] && T[i][j][1] === T[j][i][1]));
    let associative = true;
    outer: for (const a of some) for (const b of some) for (const c of some) {
      if (!same(mul(mul(a, b), c), mul(a, mul(b, c)))) {
        associative = false;
        break outer;
      }
    }
    const xs = [...some, ...pairs.slice(0, 240)];
    const alternative = associative || xs.every((x) => some.every((y) => same(mul(mul(x, x), y), mul(x, mul(x, y))) && same(mul(y, mul(x, x)), mul(mul(y, x), x))));
    const flexible = alternative || xs.every((x) => some.every((y) => same(mul(mul(x, y), x), mul(x, mul(y, x)))));
    let zeroDivisors = null;
    if (!alternative) {
      search: for (const a of pairs) for (const b of pairs) {
        if (mul(a, b).size === 0) {
          const lift = (m) => this.of(Array.from({ length: dim }, (_, k) => m.get(k) ?? 0));
          zeroDivisors = [lift(a), lift(b)];
          break search;
        }
      }
    }
    // N(xy) = N(x)N(y): true exactly when there are no zero divisors among these algebras
    const composition = alternative;
    return { name: this.name, dimension: dim, commutative, associative, alternative, flexible, composition, zeroDivisors };
  }
  toString() {
    return `${this.name} over ${this.field.name ?? "a field"}`;
  }
}

/** An element of a Cayley–Dickson algebra: `coeffs` are its coordinates on the basis 1, e1, e2, … */
export class Hyper {
  constructor(algebra, coeffs) {
    this.algebra = algebra;
    this.coeffs = Object.freeze(coeffs);
    Object.freeze(this);
  }
  #lift(o) {
    return o instanceof Hyper && o.algebra === this.algebra ? o : this.algebra.of(o);
  }
  add(o) {
    o = this.#lift(o);
    const K = this.algebra.field;
    return new Hyper(this.algebra, this.coeffs.map((c, i) => K.add(c, o.coeffs[i])));
  }
  sub(o) {
    o = this.#lift(o);
    const K = this.algebra.field;
    return new Hyper(this.algebra, this.coeffs.map((c, i) => K.sub(c, o.coeffs[i])));
  }
  neg() {
    const K = this.algebra.field;
    return new Hyper(this.algebra, this.coeffs.map((c) => K.neg(c)));
  }
  /** The product, in this order: these numbers do not commute from the quaternions on. */
  mul(o) {
    o = this.#lift(o);
    const { field: K, table: T, dim } = this.algebra;
    const out = Array(dim).fill(K.zero);
    for (let i = 0; i < dim; i++) {
      const a = this.coeffs[i];
      if (K.sign(a) === 0) continue;
      for (let j = 0; j < dim; j++) {
        const b = o.coeffs[j];
        if (K.sign(b) === 0) continue;
        const [k, s] = T[i][j];
        const p = K.mul(a, b);
        out[k] = s > 0 ? K.add(out[k], p) : K.sub(out[k], p);
      }
    }
    return new Hyper(this.algebra, out);
  }
  /** Multiply by a scalar of the field. */
  scale(c) {
    const K = this.algebra.field;
    c = (K.from ?? ((x) => x))(c);
    return new Hyper(this.algebra, this.coeffs.map((x) => K.mul(x, c)));
  }
  /** The conjugate: the same real part, every other coordinate negated. */
  conj() {
    const K = this.algebra.field;
    return new Hyper(this.algebra, this.coeffs.map((c, i) => (i === 0 ? c : K.neg(c))));
  }
  /** N(x) = x·x* = the sum of the squares of the coordinates, a scalar. */
  norm() {
    const K = this.algebra.field;
    return this.coeffs.reduce((s, c) => K.add(s, K.mul(c, c)), K.zero);
  }
  /** The real part (a scalar) and the rest. */
  get real() {
    return this.coeffs[0];
  }
  get imaginary() {
    const K = this.algebra.field;
    return new Hyper(this.algebra, this.coeffs.map((c, i) => (i === 0 ? K.zero : c)));
  }
  /** x⁻¹ = conj(x)/N(x). It satisfies x·x⁻¹ = x⁻¹·x = 1 at every level. */
  inv() {
    const K = this.algebra.field;
    const n = this.norm();
    if (K.sign(n) === 0) refuse("division-by-zero", "zero has no inverse");
    return this.conj().scale(K.div(K.one, n));
  }
  /** The x with x·b = this, that is this·b⁻¹ (unique up to the octonions; checked beyond). */
  rightDivide(b) {
    b = this.#lift(b);
    const x = this.mul(b.inv());
    if (this.algebra.dim > 8 && !x.mul(b).equals(this)) refuse("no-quotient", "this division has no answer of the form a·b⁻¹: the algebra is not alternative");
    return x;
  }
  /** The x with b·x = this, that is b⁻¹·this. */
  leftDivide(b) {
    b = this.#lift(b);
    const x = b.inv().mul(this);
    if (this.algebra.dim > 8 && !b.mul(x).equals(this)) refuse("no-quotient", "this division has no answer of the form b⁻¹·a: the algebra is not alternative");
    return x;
  }
  /** A whole-number power (well defined at every level: these algebras are power-associative). */
  pow(n) {
    if (!Number.isInteger(n)) throw new RangeError("a whole-number power is needed");
    if (n < 0) return this.inv().pow(-n);
    let result = this.algebra.one;
    let base = this;
    while (n > 0) {
      if (n & 1) result = result.mul(base);
      n >>= 1;
      if (n) base = base.mul(base);
    }
    return result;
  }
  /** ab − ba. */
  commutator(o) {
    o = this.#lift(o);
    return this.mul(o).sub(o.mul(this));
  }
  /** (ab)c − a(bc). */
  associator(b, c) {
    b = this.#lift(b);
    c = this.#lift(c);
    return this.mul(b).mul(c).sub(this.mul(b.mul(c)));
  }
  isZero() {
    const K = this.algebra.field;
    return this.coeffs.every((c) => K.sign(c) === 0);
  }
  equals(o) {
    return this.sub(this.#lift(o)).isZero();
  }
  toString() {
    const { field: K, labels } = this.algebra;
    const f = K.format ?? String;
    let out = "";
    this.coeffs.forEach((c, i) => {
      const s = K.sign(c);
      if (s === 0) return;
      const negative = s < 0 || (s !== 1 && /^-/.test(f(c)));
      const mag = f(negative ? K.neg(c) : c);
      const isOne = mag === "1";
      const wrapped = /[+\- ]/.test(mag) ? `(${mag})` : mag;
      const body = i === 0 ? wrapped : isOne ? labels[i] : `${wrapped}${labels[i].length > 1 ? " " : ""}${labels[i]}`;
      out += out ? ` ${negative ? "-" : "+"} ${body}` : `${negative ? "-" : ""}${body}`;
    });
    return out || "0";
  }
  toJSON() {
    return this.coeffs.map((c) => (this.algebra.field.format ?? String)(c));
  }
}

const cache = new Map();
/** The algebra of dimension 2^level over a field (the rationals by default). Levels 0 to 4: reals, complex, quaternions, octonions, sedenions. */
export function cayleyDickson(level, field = QQ) {
  if (field !== QQ) return new CayleyDickson(level, field);
  if (!cache.has(level)) cache.set(level, new CayleyDickson(level, QQ));
  return cache.get(level);
}
/** The quaternions, octonions and sedenions over the rationals. */
export const quaternions = cayleyDickson(2);
export const octonions = cayleyDickson(3);
export const sedenions = cayleyDickson(4);

/** The quaternion w + xi + yj + zk, or one read from text: quaternion("1 + 2i - k"). */
export const quaternion = (w = 0, x = 0, y = 0, z = 0) => (typeof w === "string" && /[ijk]/.test(w) ? quaternions.parse(w) : quaternions.of(w, x, y, z));
/** The octonion with the given eight coordinates, or one read from text: octonion("e1 + 2 e7"). */
export const octonion = (...c) => (c.length === 1 && typeof c[0] === "string" && /e\d/.test(c[0]) ? octonions.parse(c[0]) : octonions.of(...c));

// ---- quaternions as rotations ----------------------------------------------------------------

const big = (x) => (typeof x === "bigint" ? x : BigInt(x));
const absB = (x) => (x < 0n ? -x : x);
const asQuat = (q) => (q instanceof Hyper && q.algebra.dim === 4 ? q : Array.isArray(q) ? quaternions.of(q) : typeof q === "string" ? quaternions.parse(q) : quaternions.of(q));
const vec = (v) => v.map(Q);

/**
 * The rotation v → q·v·q⁻¹ as a 3×3 matrix of Rationals (rows of an array).
 * q need not have norm 1: every nonzero rational quaternion gives a rotation
 * matrix with rational entries, with no square roots anywhere.
 */
export function rotationMatrix(q) {
  const [w, x, y, z] = asQuat(q).coeffs;
  const n = asQuat(q).norm();
  if (n.numerator === 0n) refuse("division-by-zero", "the zero quaternion is not a rotation");
  const m = (a, b) => a.multiply(b);
  const two = new Rational(2n);
  const d = (a) => a.divide(n);
  return [
    [d(m(w, w).add(m(x, x)).subtract(m(y, y)).subtract(m(z, z))), d(two.multiply(m(x, y).subtract(m(w, z)))), d(two.multiply(m(x, z).add(m(w, y))))],
    [d(two.multiply(m(x, y).add(m(w, z)))), d(m(w, w).subtract(m(x, x)).add(m(y, y)).subtract(m(z, z))), d(two.multiply(m(y, z).subtract(m(w, x))))],
    [d(two.multiply(m(x, z).subtract(m(w, y)))), d(two.multiply(m(y, z).add(m(w, x)))), d(m(w, w).subtract(m(x, x)).subtract(m(y, y)).add(m(z, z)))],
  ];
}

/** Rotate the vector [x, y, z] by q: the vector part of q·v·q⁻¹. */
export function rotate(q, v) {
  q = asQuat(q);
  const [a, b, c] = vec(v);
  const r = q.mul(quaternions.of(ZERO, a, b, c)).mul(q.inv());
  return r.coeffs.slice(1);
}

/**
 * The quaternion 1 + t·(axis) for the rotation about `axis` by the angle θ
 * with tan(θ/2) = t·|axis|. For a unit axis, t is the tangent of the half
 * angle, so t = 1 is a quarter turn.
 */
export function fromAxis(axis, t = 1) {
  const [a, b, c] = vec(axis);
  t = Q(t);
  return quaternions.of(ONE, a.multiply(t), b.multiply(t), c.multiply(t));
}

/** What rotation q is: { axis, cos, tanHalfSquared } with cos θ and tan²(θ/2) exact. The axis is not normalised. */
export function rotationOf(q) {
  q = asQuat(q);
  const [w, x, y, z] = q.coeffs;
  const v2 = x.multiply(x).add(y.multiply(y)).add(z.multiply(z));
  const w2 = w.multiply(w);
  return { axis: [x, y, z], cos: w2.subtract(v2).divide(w2.add(v2)), tanHalfSquared: w2.numerator === 0n ? null : v2.divide(w2) };
}

const gcdB = (a, b) => {
  a = absB(a);
  b = absB(b);
  while (b) [a, b] = [b, a % b];
  return a;
};
/**
 * The same rotation written with whole numbers that share no factor, the first
 * nonzero one positive: the simplest quaternion on the line through q.
 */
export function primitive(q) {
  q = asQuat(q);
  if (q.isZero()) return q;
  const den = q.coeffs.reduce((l, c) => (l * c.denominator) / gcdB(l, c.denominator), 1n);
  const ints = q.coeffs.map((c) => (c.numerator * den) / c.denominator);
  const g = ints.reduce(gcdB, 0n);
  const lead = ints.find((x) => x !== 0n);
  return quaternions.of(ints.map((x) => new Rational((lead < 0n ? -x : x) / g)));
}

/**
 * A quaternion for a rotation matrix with rational entries (orthogonal,
 * determinant 1), found without square roots. The answer is determined up to
 * a scalar factor (the primitive one is returned), and rotationMatrix gives the matrix back.
 */
export function fromRotationMatrix(M) {
  const m = M.map((row) => vec(row));
  const dot = (a, b) => a.reduce((s, x, i) => s.add(x.multiply(b[i])), ZERO);
  for (let i = 0; i < 3; i++) for (let j = i; j < 3; j++) if (!dot(m[i], m[j]).equals(i === j ? ONE : ZERO)) refuse("not-a-rotation", "the matrix is not orthogonal");
  const tr = m[0][0].add(m[1][1]).add(m[2][2]);
  // four candidates, proportional to the same quaternion; take one that is not zero
  const cands = [
    [ONE.add(tr), m[2][1].subtract(m[1][2]), m[0][2].subtract(m[2][0]), m[1][0].subtract(m[0][1])],
    [m[2][1].subtract(m[1][2]), ONE.add(m[0][0]).subtract(m[1][1]).subtract(m[2][2]), m[0][1].add(m[1][0]), m[0][2].add(m[2][0])],
    [m[0][2].subtract(m[2][0]), m[0][1].add(m[1][0]), ONE.subtract(m[0][0]).add(m[1][1]).subtract(m[2][2]), m[1][2].add(m[2][1])],
    [m[1][0].subtract(m[0][1]), m[0][2].add(m[2][0]), m[1][2].add(m[2][1]), ONE.subtract(m[0][0]).subtract(m[1][1]).add(m[2][2])],
  ];
  const c = cands.find((k) => k.some((x) => x.numerator !== 0n));
  const q = primitive(quaternions.of(c));
  const back = rotationMatrix(q);
  if (!back.every((row, i) => row.every((x, j) => x.equals(m[i][j])))) refuse("not-a-rotation", "the matrix is orthogonal but has determinant −1 (a reflection)");
  return q;
}

/**
 * A rotation part-way from q0 to q1: (1 − t)·q0 + t·q1, with q1's sign chosen
 * for the shorter way round. Every value is an exact rational rotation and
 * the path is the same great arc that spherical interpolation follows (at a
 * different speed).
 */
export function interpolate(q0, q1, t) {
  q0 = asQuat(q0);
  q1 = asQuat(q1);
  t = Q(t);
  // compare directions after scaling: the sign of the dot product decides the shorter arc
  const dot = q0.coeffs.reduce((s, c, i) => s.add(c.multiply(q1.coeffs[i])), ZERO);
  if (dot.numerator < 0n) q1 = q1.neg();
  return q0.scale(ONE.subtract(t)).add(q1.scale(t));
}

// ---- integer quaternions ---------------------------------------------------------------------
// A Hurwitz quaternion has coordinates that are all whole numbers or all
// halves of odd numbers. It is held as four BigInts [a, b, c, d] meaning
// (a + bi + cj + dk)/2, all of the same parity. Lipschitz quaternions are the
// ones with whole-number coordinates (a, b, c, d all even).

const sameParity = (v) => v.every((x) => ((x % 2n) + 2n) % 2n === ((v[0] % 2n) + 2n) % 2n);

export class Hurwitz {
  /** From doubled coordinates: new Hurwitz([1n, 1n, 1n, 1n]) is (1 + i + j + k)/2. */
  constructor(doubled) {
    const v = doubled.map(big);
    if (v.length !== 4 || !sameParity(v)) throw new RangeError("a Hurwitz quaternion has coordinates that are all whole numbers or all halves of odd numbers");
    this.d = Object.freeze(v);
    Object.freeze(this);
  }
  /** From a Hurwitz, a rational quaternion, or four coordinates that are whole numbers or halves ("1/2"). */
  static from(x, ...rest) {
    if (x instanceof Hurwitz) return x;
    const cs = x instanceof Hyper ? x.coeffs : Array.isArray(x) ? x : [x, ...rest];
    const q = Array.from({ length: 4 }, (_, i) => (cs[i] === undefined ? ZERO : Q(cs[i])).multiply(new Rational(2n)));
    if (q.some((c) => c.denominator !== 1n)) throw new RangeError("a Hurwitz quaternion has coordinates that are whole numbers or halves");
    return new Hurwitz(q.map((c) => c.numerator));
  }
  /** Whole-number coordinates. */
  get isLipschitz() {
    return this.d[0] % 2n === 0n;
  }
  /** The norm a² + b² + c² + d², a whole number. */
  norm() {
    return this.d.reduce((s, x) => s + x * x, 0n) / 4n;
  }
  conj() {
    const [a, b, c, d] = this.d;
    return new Hurwitz([a, -b, -c, -d]);
  }
  neg() {
    return new Hurwitz(this.d.map((x) => -x));
  }
  add(o) {
    o = Hurwitz.from(o);
    return new Hurwitz(this.d.map((x, i) => x + o.d[i]));
  }
  sub(o) {
    o = Hurwitz.from(o);
    return new Hurwitz(this.d.map((x, i) => x - o.d[i]));
  }
  mul(o) {
    o = Hurwitz.from(o);
    const [a, b, c, d] = this.d;
    const [e, f, g, h] = o.d;
    return new Hurwitz([a * e - b * f - c * g - d * h, a * f + b * e + c * h - d * g, a * g - b * h + c * e + d * f, a * h + b * g - c * f + d * e].map((x) => x / 2n));
  }
  /**
   * Division with remainder on the right: { quotient, remainder } with
   * this = quotient·b + remainder and N(remainder) < N(b). It always exists
   * for Hurwitz quaternions, which is what makes the Euclidean algorithm work.
   */
  divmod(b) {
    b = Hurwitz.from(b);
    const n = b.norm();
    if (n === 0n) refuse("division-by-zero", "division by zero");
    // this·b* / N(b), then the nearest Hurwitz quaternion: the nearer of the nearest Lipschitz point and the nearest half-odd point
    const t = this.mul(b.conj()).d; // doubled coordinates of this·b*, to be divided by n
    const roundTo = (parity) => t.map((x) => {
      // nearest integer m ≡ parity (mod 2) to x/n
      const f = floorDiv(x - parity * n, 2n * n);
      const lo = 2n * f + parity;
      return absB(x - lo * n) <= absB(x - (lo + 2n) * n) ? lo : lo + 2n;
    });
    const dist = (m) => m.reduce((s, v, i) => s + (t[i] - v * n) ** 2n, 0n);
    const even = roundTo(0n);
    const odd = roundTo(1n);
    const q = new Hurwitz(dist(even) <= dist(odd) ? even : odd);
    return { quotient: q, remainder: this.sub(q.mul(b)) };
  }
  isZero() {
    return this.d.every((x) => x === 0n);
  }
  equals(o) {
    o = Hurwitz.from(o);
    return this.d.every((x, i) => x === o.d[i]);
  }
  /** As a quaternion over the rationals. */
  toQuaternion() {
    return quaternions.of(this.d.map((x) => new Rational(x, 2n)));
  }
  /** The whole-number coordinates [a, b, c, d] of a Lipschitz quaternion. */
  toIntegers() {
    if (!this.isLipschitz) throw new RangeError("this Hurwitz quaternion has half-integer coordinates");
    return this.d.map((x) => x / 2n);
  }
  toString() {
    return this.toQuaternion().toString();
  }
}
const floorDiv = (a, b) => {
  const q = a / b;
  return a % b !== 0n && (a < 0n) !== (b < 0n) ? q - 1n : q;
};

/** hurwitz(1, 2, 0, -1) or hurwitz("1/2", "1/2", "1/2", "1/2"). */
export const hurwitz = (...c) => Hurwitz.from(...c);

/** The 24 units: ±1, ±i, ±j, ±k and (±1 ± i ± j ± k)/2. */
export function hurwitzUnits() {
  const out = [];
  for (let i = 0; i < 4; i++) for (const s of [2n, -2n]) out.push(new Hurwitz([0n, 0n, 0n, 0n].map((_, k) => (k === i ? s : 0n))));
  for (let m = 0; m < 16; m++) out.push(new Hurwitz([0, 1, 2, 3].map((k) => ((m >> k) & 1 ? -1n : 1n))));
  return out;
}

/** A right greatest common divisor: g with a = a′·g and b = b′·g, by the Euclidean algorithm. Determined up to a unit on the left. */
export function hurwitzGcd(a, b) {
  a = Hurwitz.from(a);
  b = Hurwitz.from(b);
  while (!b.isZero()) [a, b] = [b, a.divmod(b).remainder];
  return a;
}

const modpow = (b, e, m) => {
  let r = 1n;
  b %= m;
  while (e > 0n) {
    if (e & 1n) r = (r * b) % m;
    b = (b * b) % m;
    e >>= 1n;
  }
  return r;
};
/** Miller–Rabin with the first twelve primes as bases: exact below 3.3·10²⁴. */
function isPrime(n) {
  if (n < 2n) return false;
  const small = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n];
  for (const p of small) if (n % p === 0n) return n === p;
  let d = n - 1n;
  let s = 0;
  while (d % 2n === 0n) {
    d /= 2n;
    s++;
  }
  next: for (const a of small) {
    let x = modpow(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    for (let i = 1; i < s; i++) {
      x = (x * x) % n;
      if (x === n - 1n) continue next;
    }
    return false;
  }
  return true;
}
/** A square root of a modulo an odd prime p (Tonelli–Shanks), or null. */
function sqrtMod(a, p) {
  a = ((a % p) + p) % p;
  if (a === 0n) return 0n;
  if (modpow(a, (p - 1n) / 2n, p) !== 1n) return null;
  if (p % 4n === 3n) return modpow(a, (p + 1n) / 4n, p);
  let q = p - 1n;
  let s = 0n;
  while (q % 2n === 0n) {
    q /= 2n;
    s++;
  }
  let z = 2n;
  while (modpow(z, (p - 1n) / 2n, p) !== p - 1n) z++;
  let m = s;
  let c = modpow(z, q, p);
  let t = modpow(a, q, p);
  let r = modpow(a, (q + 1n) / 2n, p);
  while (t !== 1n) {
    let i = 0n;
    let u = t;
    while (u !== 1n) {
      u = (u * u) % p;
      i++;
    }
    const b = modpow(c, 1n << (m - i - 1n), p);
    m = i;
    c = (b * b) % p;
    t = (t * c) % p;
    r = (r * b) % p;
  }
  return r;
}
function factorSmall(n) {
  const out = [];
  let p = 2n;
  while (p * p <= n && p < 1000000n) {
    while (n % p === 0n) {
      out.push(p);
      n /= p;
    }
    p += p === 2n ? 1n : 2n;
  }
  if (n > 1n) {
    if (!isPrime(n)) refuse("too-large", `cannot factor ${n} here: it has no prime factor below a million and is not prime`);
    out.push(n);
  }
  return out;
}

/** A Lipschitz quaternion of norm p for a prime p: the heart of Lagrange's theorem. */
export function quaternionOfNorm(p) {
  p = big(p);
  if (!isPrime(p)) throw new RangeError(`${p} is not prime`);
  if (p === 2n) return new Hurwitz([2n, 2n, 0n, 0n]);
  // x² + y² + 1 ≡ 0 (mod p) always has a solution; then p divides N(x + yi + j)
  let x = 0n;
  let y = null;
  for (;; x++) {
    y = sqrtMod(-1n - x * x, p);
    if (y !== null) break;
  }
  // the right gcd of x + yi + j and p has norm exactly p
  let g = hurwitzGcd(new Hurwitz([2n * x, 2n * y, 2n, 0n]), new Hurwitz([2n * p, 0n, 0n, 0n]));
  if (!g.isLipschitz) g = hurwitzUnits().map((u) => u.mul(g)).find((h) => h.isLipschitz);
  return g;
}

/**
 * Four whole numbers whose squares add to n (Lagrange's theorem), found the
 * way the theorem is proven: a quaternion of norm p for each prime factor,
 * multiplied together, since the norm is multiplicative. Returned largest
 * first, as BigInts.
 */
export function fourSquares(n) {
  n = big(n);
  if (n < 0n) throw new RangeError("a sum of squares is not negative");
  if (n === 0n) return [0n, 0n, 0n, 0n];
  let q = new Hurwitz([2n, 0n, 0n, 0n]);
  for (const p of factorSmall(n)) q = q.mul(quaternionOfNorm(p));
  return q.toIntegers().map(absB).sort((a, b) => (a < b ? 1 : a > b ? -1 : 0));
}

/** The number of ways to write n as an ordered sum of four squares of integers (Jacobi): 8 times the sum of the divisors of n not divisible by 4. */
export function fourSquaresCount(n) {
  n = big(n);
  if (n < 0n) return 0n;
  if (n === 0n) return 1n;
  const ps = factorSmall(n);
  const exps = new Map();
  for (const p of ps) exps.set(p, (exps.get(p) ?? 0) + 1);
  let divisors = [1n];
  for (const [p, e] of exps) {
    const next = [];
    for (const d of divisors) for (let k = 0, pk = 1n; k <= e; k++, pk *= p) next.push(d * pk);
    divisors = next;
  }
  return 8n * divisors.filter((d) => d % 4n !== 0n).reduce((s, d) => s + d, 0n);
}
