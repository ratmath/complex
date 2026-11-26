// Möbius transformations, generalized circles, the Riemann sphere and
// conformal maps of standard regions.
//
// A Möbius map z ↦ (az + b)/(cz + d) with complex rational a, b, c, d is
// exact: it sends complex rationals to complex rationals and the point at
// infinity where it should. Circles and lines (generalized circles) are
// kept as Hermitian forms with rational data, and a Möbius map sends one to
// another exactly, so the image of a grid of lines is a family of exact
// circles. Stereographic projection to the Riemann sphere is rational too.
//
// Fixed points need a square root, so they are exact when the root is a
// complex rational and certified rectangles otherwise.
import { Rational } from "@ratmath/core";
import { Complex, ZERO, ONE, TWO, refuse } from "./complex.js";
import { Rect } from "./enclosure.js";
import { sqrt as sqrtB } from "./functions.js";

/** The point at infinity of the Riemann sphere. */
export const infinity = Object.freeze({
  infinity: true,
  toString: () => "∞",
  toJSON: () => "infinity",
});
export const isInfinity = (z) => z === infinity || z === "infinity" || z === Infinity;
const C = (z) => (isInfinity(z) ? infinity : Complex.from(z));
const czero = Complex.ZERO;
const cone = Complex.ONE;
const isZ = (z) => z.isZero();
const eqPt = (z, w) => (z === infinity || w === infinity ? z === w : z.equals(w));

/** A Möbius transformation z ↦ (az + b)/(cz + d), ad − bc ≠ 0. */
export class Mobius {
  constructor(a, b, c, d) {
    [this.a, this.b, this.c, this.d] = [a, b, c, d].map((x) => Complex.from(x));
    if (isZ(this.det)) refuse("degenerate", "ad − bc = 0, so this is not a Möbius map (it is constant)");
    Object.freeze(this);
  }
  get det() {
    return this.a.mul(this.d).sub(this.b.mul(this.c));
  }
  /** The matrix [[a, b], [c, d]]. */
  get matrix() {
    return [[this.a, this.b], [this.c, this.d]];
  }
  /** Where z goes: a Complex, or infinity. */
  apply(z) {
    const { a, b, c, d } = this;
    if (isInfinity(z)) return isZ(c) ? infinity : a.div(c);
    const w = Complex.from(z);
    const den = c.mul(w).add(d);
    if (isZ(den)) return infinity;
    return a.mul(w).add(b).div(den);
  }
  /** this after other: z ↦ this(other(z)). */
  compose(o) {
    const [[a, b], [c, d]] = this.matrix;
    const [[e, f], [g, h]] = o.matrix;
    return new Mobius(a.mul(e).add(b.mul(g)), a.mul(f).add(b.mul(h)), c.mul(e).add(d.mul(g)), c.mul(f).add(d.mul(h)));
  }
  inverse() {
    return new Mobius(this.d, this.b.neg(), this.c.neg(), this.a);
  }
  /** Equal as maps: the matrices are proportional. */
  equals(o) {
    const P = this.matrix.flat();
    const R = o.matrix.flat();
    // P ∝ R ⇔ every 2 × 2 minor of [P; R] vanishes
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) if (!P[i].mul(R[j]).sub(P[j].mul(R[i])).isZero()) return false;
    return true;
  }
  /** (a + d)²/(ad − bc): unchanged when the matrix is scaled, and it decides the kind of map. */
  get traceSquared() {
    const t = this.a.add(this.d);
    return t.mul(t).div(this.det);
  }
  get isIdentity() {
    return this.equals(new Mobius(1, 0, 0, 1));
  }
  /**
   * The kind of map, from τ = (a + d)²/(ad − bc): identity; parabolic (τ = 4,
   * one fixed point); elliptic (τ real in [0, 4): a rotation about two fixed
   * points); hyperbolic (τ real > 4: a pure stretch between them); loxodromic
   * (anything else: a stretch and a turn).
   */
  classify() {
    if (this.isIdentity) return "identity";
    const tau = this.traceSquared;
    if (!tau.im.equals(ZERO)) return "loxodromic";
    const t = tau.re;
    if (t.equals(new Rational(4n))) return "parabolic";
    if (t.compareTo(ZERO) >= 0 && t.compareTo(new Rational(4n)) < 0) return "elliptic";
    if (t.compareTo(new Rational(4n)) > 0) return "hyperbolic";
    return "loxodromic";
  }
  /**
   * The fixed points: roots of cz² + (d − a)z − b = 0 (with ∞ when c = 0).
   * { exact, points }: Complex values (or infinity) when the square root is a
   * complex rational, else Rect enclosures.
   */
  fixedPoints({ bits = 64 } = {}) {
    const { a, b, c, d } = this;
    if (this.isIdentity) return Object.freeze({ exact: true, points: "every point" });
    if (isZ(c)) {
      if (a.equals(d)) return Object.freeze({ exact: true, points: [infinity] });
      return Object.freeze({ exact: true, points: [b.div(d.sub(a)), infinity] });
    }
    const amd = a.sub(d);
    const disc = amd.mul(amd).add(Complex.from(4).mul(b).mul(c));
    const two_c = c.mul(Complex.from(2));
    if (disc.isZero()) return Object.freeze({ exact: true, points: [amd.div(two_c)] });
    const w = disc.exactSqrt();
    if (w) return Object.freeze({ exact: true, points: [amd.add(w[0]).div(two_c), amd.sub(w[0]).div(two_c)] });
    // either square root will do (the two fixed points use both signs); on the cut take one side
    const root = sqrtB(disc, { bits });
    const s = root.value ?? root.candidates[0].value;
    const A = Rect.from(amd);
    const D = Rect.from(two_c);
    return Object.freeze({ exact: false, points: [A.add(s).div(D), A.sub(s).div(D)] });
  }
  /** Is this a rotation of the Riemann sphere? (Proportional to a unitary matrix: d = ā, c = −b̄ up to a common factor.) */
  get isSphereRotation() {
    const { a, b, c, d } = this;
    // M M* ∝ I: |a|² + |b|² = |c|² + |d|² and a c̄ + b d̄ = 0
    return a.norm().add(b.norm()).equals(c.norm().add(d.norm())) && a.mul(c.conj()).add(b.mul(d.conj())).isZero();
  }
  /** The map on doubles, for pictures: [x, y] → [u, v] (null at infinity). */
  applyFloat([x, y]) {
    const [a, b, c, d] = [this.a, this.b, this.c, this.d].map((z) => z.toNumbers());
    const nr = a[0] * x - a[1] * y + b[0];
    const ni = a[0] * y + a[1] * x + b[1];
    const dr = c[0] * x - c[1] * y + d[0];
    const di = c[0] * y + c[1] * x + d[1];
    const m = dr * dr + di * di;
    if (m === 0) return null;
    return [(nr * dr + ni * di) / m, (ni * dr - nr * di) / m];
  }
  toString() {
    const p = (z) => (z.re.equals(ZERO) || z.im.equals(ZERO) ? String(z) : `(${z})`);
    const term = (k, c) => (isZ(k) ? null : k.equals(cone) ? "z" : `${p(k)}z`);
    const side = (k, c) => [term(k), isZ(c) ? null : p(c)].filter(Boolean).join(" + ") || "0";
    return `z ↦ (${side(this.a, this.b)})/(${side(this.c, this.d)})`.replace(/\+ -/g, "− ");
  }
}
export const mobius = (a, b, c, d) => new Mobius(a, b, c, d);

/** The map z ↦ (z − z₁)(z₂ − z₃)/((z − z₃)(z₂ − z₁)) sending z₁, z₂, z₃ to 0, 1, ∞ (any one may be ∞). */
function toStandard(z1, z2, z3) {
  const [p, q, r] = [z1, z2, z3].map(C);
  if (eqPt(p, q) || eqPt(q, r) || eqPt(p, r)) refuse("degenerate", "the three points must be distinct");
  if (p === infinity) return new Mobius(czero, q.sub(r), cone, r.neg());
  if (q === infinity) return new Mobius(cone, p.neg(), cone, r.neg());
  if (r === infinity) return new Mobius(cone, p.neg(), czero, q.sub(p));
  // (z − p)(q − r) / ((z − r)(q − p))
  return new Mobius(q.sub(r), p.neg().mul(q.sub(r)), q.sub(p), r.neg().mul(q.sub(p)));
}
/** The one Möbius map sending three distinct points to three distinct points. */
export function mobiusFrom(from, to) {
  if (from.length !== 3 || to.length !== 3) refuse("bad-argument", "give three points and their three images");
  return toStandard(...to).inverse().compose(toStandard(...from));
}
/**
 * The cross-ratio (z₁, z₂; z₃, z₄) = (z₁ − z₃)(z₂ − z₄)/((z₂ − z₃)(z₁ − z₄)):
 * the image of z₁ under the map sending z₂, z₃, z₄ to 1, 0, ∞. Möbius maps
 * keep it, and it is real exactly when the four points lie on one circle or line.
 */
export function crossRatio(z1, z2, z3, z4) {
  // (z₁, z₂; z₃, z₄) = S(z₁) for S sending z₃ → 0, z₂ → 1, z₄ → ∞
  return toStandard(z3, z2, z4).apply(C(z1));
}
/** Do four points lie on one circle or line? Exactly when the cross-ratio is real (or ∞). */
export function concyclic(z1, z2, z3, z4) {
  const r = crossRatio(z1, z2, z3, z4);
  return r === infinity || r.im.equals(ZERO);
}

// ---- generalized circles ------------------------------------------------------------

/**
 * A circle or a line, as A|z|² + B̄z + Bz̄ + C = 0 with A, C rational and B
 * complex rational (the Hermitian matrix [[A, B], [B̄, C]]). A ≠ 0: the
 * circle with centre −B/A and radius² (|B|² − AC)/A². A = 0: the line.
 */
export class GCircle {
  constructor(A, B, C) {
    this.A = Complex.from(A).re;
    this.B = Complex.from(B);
    this.C = Complex.from(C).re;
    if (!Complex.from(A).im.equals(ZERO) || !Complex.from(C).im.equals(ZERO)) refuse("bad-argument", "A and C are real");
    if (this.radius2 !== null && this.radius2.compareTo(ZERO) <= 0) refuse("degenerate", "this circle has no real points (radius² ≤ 0)");
    if (this.A.equals(ZERO) && this.B.isZero()) refuse("degenerate", "not a circle or a line");
    Object.freeze(this);
  }
  get isLine() {
    return this.A.equals(ZERO);
  }
  get center() {
    return this.isLine ? null : this.B.neg().scale(ONE.divide(this.A));
  }
  get radius2() {
    if (this.A.equals(ZERO)) return null;
    return this.B.norm().subtract(this.A.multiply(this.C)).divide(this.A.multiply(this.A));
  }
  /** The value of the form at z (0 on the circle; the sign tells the side). Lines contain ∞. */
  value(z) {
    if (isInfinity(z)) return this.A;
    const w = Complex.from(z);
    // A|z|² + 2 Re(B̄ z) + C
    return this.A.multiply(w.norm()).add(TWO.multiply(this.B.conj().mul(w).re)).add(this.C);
  }
  contains(z) {
    return this.value(z).equals(ZERO);
  }
  /** Same circle: proportional forms. */
  equals(o) {
    const P = [this.A, this.B.re, this.B.im, this.C];
    const R = [o.A, o.B.re, o.B.im, o.C];
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) if (!P[i].multiply(R[j]).equals(P[j].multiply(R[i]))) return false;
    return true;
  }
  /** The image under a Möbius map T: H′ = (T⁻¹)* H T⁻¹. */
  map(T) {
    const [[p, q], [r, s]] = T.inverse().matrix; // T⁻¹ = [[p, q], [r, s]]
    // z = (p w + q)/(r w + s); the form in w is the conjugate transpose sandwich
    const H = [[Complex.from(this.A), this.B], [this.B.conj(), Complex.from(this.C)]];
    const M = [[p, q], [r, s]];
    const Ms = [[p.conj(), r.conj()], [q.conj(), s.conj()]];
    const mul = (X, Y) => [[X[0][0].mul(Y[0][0]).add(X[0][1].mul(Y[1][0])), X[0][0].mul(Y[0][1]).add(X[0][1].mul(Y[1][1]))], [X[1][0].mul(Y[0][0]).add(X[1][1].mul(Y[1][0])), X[1][0].mul(Y[0][1]).add(X[1][1].mul(Y[1][1]))]];
    const K = mul(Ms, mul(H, M));
    return new GCircle(K[0][0].re, K[0][1], K[1][1].re);
  }
  toString() {
    if (this.isLine) {
      // 2 Re(B̄ z) + C = 0 is 2Bx·x + 2By·y + C = 0
      const terms = [[TWO.multiply(this.B.re), "x"], [TWO.multiply(this.B.im), "y"], [this.C, ""]].filter(([c]) => !c.equals(ZERO));
      const text = terms.map(([c, v], i) => {
        const neg = c.compareTo(ZERO) < 0;
        const m = neg ? c.negate() : c;
        const coef = v && m.equals(ONE) ? "" : String(m);
        return `${i ? (neg ? " − " : " + ") : neg ? "−" : ""}${coef}${v}`;
      }).join("");
      return `the line ${text} = 0`;
    }
    const c = this.center;
    return `the circle ${c.isZero() ? "|z|" : `|z − (${c})|`}² = ${this.radius2}`;
  }
}
/** The circle with centre c and radius² r2 (rational), exactly. */
export function circle(center, radius2) {
  const c = Complex.from(center);
  const r2 = Complex.from(radius2).re;
  // |z − c|² − r² = |z|² − c̄z − cz̄ + |c|² − r²: A = 1, B = −c, C = |c|² − r²
  return new GCircle(ONE, c.neg(), c.norm().subtract(r2));
}
/** The generalized circle through three points (a line if they are collinear or one is ∞). */
export function circleThrough(z1, z2, z3) {
  const pts = [z1, z2, z3].map(C).filter((z) => z !== infinity);
  if (pts.length < 2) refuse("degenerate", "give at least two finite points");
  if (pts.length === 2) return lineThrough(pts[0], pts[1]);
  const [a, b, c] = pts;
  const cross = b.sub(a).conj().mul(c.sub(a)).im;
  if (cross.equals(ZERO)) return lineThrough(a, b);
  // solve A|z|² + 2Re(B̄z) + C = 0 with A = 1 for the three points: linear in (Bx, By, C)
  // |z|² + 2(Bx x + By y) + C = 0
  const rows = [a, b, c].map((z) => [TWO.multiply(z.re), TWO.multiply(z.im), ONE, z.norm().negate()]);
  const sol = solve3(rows);
  return new GCircle(ONE, new Complex(sol[0], sol[1]), sol[2]);
}
/** The line through two distinct points. */
export function lineThrough(z1, z2) {
  const a = Complex.from(z1);
  const b = Complex.from(z2);
  if (a.equals(b)) refuse("degenerate", "the two points are the same");
  // a normal n = i(b − a); 2Re(B̄z) + C = 0 with B = n/2, C = −2Re(B̄a)
  const B = Complex.I.mul(b.sub(a));
  const C2 = TWO.multiply(B.conj().mul(a).re).negate();
  return new GCircle(ZERO, B, C2);
}
function solve3(rows) {
  const A = rows.map((r) => [...r]);
  for (let c = 0; c < 3; c++) {
    const p = A.findIndex((r, i) => i >= c && !r[c].equals(ZERO));
    [A[c], A[p]] = [A[p], A[c]];
    const pv = A[c][c];
    A[c] = A[c].map((x) => x.divide(pv));
    for (let i = 0; i < 3; i++) if (i !== c && !A[i][c].equals(ZERO)) {
      const f = A[i][c];
      A[i] = A[i].map((x, j) => x.subtract(f.multiply(A[c][j])));
    }
  }
  return A.map((r) => r[3]);
}

// ---- the Riemann sphere ---------------------------------------------------------------

/** Stereographic projection to the unit sphere from the north pole: z ↦ (2x, 2y, |z|² − 1)/(|z|² + 1); ∞ ↦ (0, 0, 1). Exact. */
export function toSphere(z) {
  if (isInfinity(z)) return [ZERO, ZERO, ONE];
  const w = Complex.from(z);
  const n = w.norm();
  const den = n.add(ONE);
  return [TWO.multiply(w.re).divide(den), TWO.multiply(w.im).divide(den), n.subtract(ONE).divide(den)];
}
/** The inverse: (X, Y, Z) ↦ (X + iY)/(1 − Z); the north pole goes to ∞. The point must be on the sphere. */
export function fromSphere([X, Y, Z]) {
  const [x, y, z] = [X, Y, Z].map((t) => Complex.from(t).re);
  if (!x.multiply(x).add(y.multiply(y)).add(z.multiply(z)).equals(ONE)) refuse("bad-argument", "the point is not on the unit sphere");
  if (z.equals(ONE)) return infinity;
  return new Complex(x, y).scale(ONE.divide(ONE.subtract(z)));
}
/** The chordal distance between two points of the Riemann sphere, squared: 4|z − w|²/((1 + |z|²)(1 + |w|²)). */
export function chordal2(z, w) {
  const [P, R] = [z, w].map(toSphere);
  return P.reduce((s, x, i) => s.add(x.subtract(R[i]).multiply(x.subtract(R[i]))), ZERO);
}

// ---- conformal maps of standard regions -------------------------------------------------

/** z ↦ (z − i)/(z + i): the upper half-plane onto the unit disc (the Cayley map). */
export const cayley = () => new Mobius(1, Complex.I.neg(), 1, Complex.I);
/**
 * The disc automorphism z ↦ u(z − a)/(1 − āz), |a| < 1, |u| = 1 (u exact,
 * for instance (3 + 4i)/5). It sends a to 0. Every conformal self-map of the
 * disc has this form.
 */
export function discAutomorphism(a, u = 1) {
  const A = Complex.from(a);
  const U = Complex.from(u);
  if (A.norm().compareTo(ONE) >= 0) refuse("bad-argument", "a must be inside the unit disc");
  if (!U.norm().equals(ONE)) refuse("bad-argument", "u must have modulus 1 exactly");
  return new Mobius(U, U.mul(A).neg(), A.conj().neg(), cone);
}
/** The map of the upper half-plane onto the disc sending z₀ to 0: z ↦ (z − z₀)/(z − z̄₀). */
export function halfPlaneToDisc(z0) {
  const z = Complex.from(z0);
  if (z.im.compareTo(ZERO) <= 0) refuse("bad-argument", "z₀ must be in the upper half-plane");
  return new Mobius(cone, z.neg(), cone, z.conj().neg());
}
/** The Joukowski map z ↦ (z + 1/z)/2: the outside of the unit disc onto the plane cut along [−1, 1]; circles |z| = r go to ellipses. Exact on complex rationals. */
export function joukowski(z) {
  const w = Complex.from(z);
  if (w.isZero()) return infinity;
  return w.add(w.inv()).scale(new Rational(1n, 2n));
}

/** Membership in standard regions, exactly. */
export const regions = Object.freeze({
  disc: (z) => !isInfinity(z) && Complex.from(z).norm().compareTo(ONE) < 0,
  upperHalfPlane: (z) => !isInfinity(z) && Complex.from(z).im.compareTo(ZERO) > 0,
  rightHalfPlane: (z) => !isInfinity(z) && Complex.from(z).re.compareTo(ZERO) > 0,
  firstQuadrant: (z) => !isInfinity(z) && Complex.from(z).re.compareTo(ZERO) > 0 && Complex.from(z).im.compareTo(ZERO) > 0,
});

// ---- graphics hooks ----------------------------------------------------------------------

/**
 * The image of a grid of horizontal and vertical lines under a Möbius map:
 * each line goes to an exact generalized circle. The result is data for a
 * picture: { kind: "exact", lines: [{ from: "re = x" | "im = y", image: GCircle }] }.
 */
export function gridImage(T, { re = [-2, 2], im = [-2, 2], step = "1/2" } = {}) {
  const st = Complex.from(step).re;
  const out = [];
  const range = ([lo, hi]) => {
    const vals = [];
    for (let x = Complex.from(lo).re; x.compareTo(Complex.from(hi).re) <= 0; x = x.add(st)) vals.push(x);
    if (vals.length > 400) refuse("too-large", "too many grid lines");
    return vals;
  };
  for (const x of range(re)) out.push({ from: `re = ${x}`, image: lineThrough(new Complex(x, ZERO), new Complex(x, ONE)).map(T) });
  for (const y of range(im)) out.push({ from: `im = ${y}`, image: lineThrough(new Complex(ZERO, y), new Complex(ONE, y)).map(T) });
  return Object.freeze({ kind: "exact", lines: out });
}
/**
 * Values of a map on a grid of doubles, for domain colouring or a picture of
 * the image of a grid: f is a Mobius, or any function [x, y] → [u, v]
 * (null where undefined). Marked estimated: these are floating-point values.
 */
export function sampleGrid(f, { re = [-2, 2], im = [-2, 2], width = 64, height = 64 } = {}) {
  const g = f instanceof Mobius ? (p) => f.applyFloat(p) : f;
  const [x0, x1] = re.map(Number);
  const [y0, y1] = im.map(Number);
  const values = new Float64Array(2 * width * height);
  for (let j = 0; j < height; j++)
    for (let i = 0; i < width; i++) {
      const p = [x0 + ((i + 0.5) * (x1 - x0)) / width, y0 + ((j + 0.5) * (y1 - y0)) / height];
      const v = g(p);
      const k = 2 * (j * width + i);
      values[k] = v ? v[0] : NaN;
      values[k + 1] = v ? v[1] : NaN;
    }
  return Object.freeze({ kind: "estimated", width, height, re: [x0, x1], im: [y0, y1], values });
}

