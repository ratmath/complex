// Contours, proven contour integrals, the argument principle, and zeros of
// analytic functions with their multiplicities.
//
// Everything here stands on Taylor models (taylor.js). A model on a disc
// proves the function analytic there and encloses it, so:
//
//   • the integral over a piece of contour inside the disc is the difference
//     of an antiderivative at the two ends, whatever the piece's shape;
//   • when the model's range misses 0, the function has no zero in the disc
//     and its argument changes by less than half a turn along the piece, so
//     the change is the principal argument of f(end)/f(start);
//   • a Newton step on a disc that lands inside the disc proves exactly one
//     simple zero there.
//
// Pieces are subdivided until the statements can be proved. What cannot be
// proved is reported, never guessed.

import { Rational } from "@ratmath/core";
import { constants } from "@ratmath/analytic";
import { Complex, ComplexError, Q, ZERO, ONE, TWO, refuse } from "./complex.js";
import { Rect, Disc, I, iv, sqrtUp, roundUp } from "./enclosure.js";
import { arg } from "./functions.js";
import { CTM, asFunction } from "./taylor.js";

const HALF = new Rational(1n, 2n);
const PAD = new Rational(65n, 64n);
const SHARP = new Rational(1n, 1n << 16n);
const cross = (u, v) => u.re.multiply(v.im).subtract(u.im.multiply(v.re));
const dot = (u, v) => u.re.multiply(v.re).add(u.im.multiply(v.im));
const refusal = (e) => {
  if (e instanceof ComplexError) return e;
  throw e;
};

// ---- pieces ---------------------------------------------------------------------------------

/** A straight piece, or (with a centre) a piece of circle of at most a quarter turn. */
const line = (a, b) => ({ a, b, center: null, ccw: true });

/** A rational point on the circle about c through a, about `fraction` of the way round to b. */
function arcPoint(a, b, c, ccw, fraction = 0.5) {
  const v = b.sub(c).div(a.sub(c));
  const [x, y] = v.toNumbers();
  let th = Math.atan2(y, x);
  if (ccw && th <= 0) th += 2 * Math.PI;
  if (!ccw && th >= 0) th -= 2 * Math.PI;
  // e/ē is a rational point of the unit circle at twice the angle of e, whatever the size of e
  const half = (th * fraction) / 2;
  let e;
  if (Math.abs(half) < 1e-3) {
    // a small angle: 1 + i·tan, with the tangent kept to 26 significant bits
    const t = Math.tan(half);
    const k = 26 - Math.floor(Math.log2(Math.abs(t)));
    e = new Complex(ONE, new Rational(BigInt(Math.round(t * 2 ** k)), 1n << BigInt(k)));
  } else {
    const S = 2 ** 26;
    e = new Complex(new Rational(BigInt(Math.round(Math.cos(half) * S)), BigInt(S)), new Rational(BigInt(Math.round(Math.sin(half) * S)), BigInt(S)));
  }
  return c.add(a.sub(c).mul(e.div(e.conj())));
}
const isQuarter = ({ a, b, center, ccw }) => {
  const [u, v] = [a.sub(center), b.sub(center)];
  const s = cross(u, v).numerator;
  return dot(u, v).numerator >= 0n && (ccw ? s > 0n : s < 0n);
};
function arcPieces(a, b, center, ccw, depth = 0) {
  const piece = { a, b, center, ccw };
  if (depth > 6) refuse("bad-argument", "could not cut the arc into quarter turns");
  if (isQuarter(piece)) return [piece];
  const m = arcPoint(a, b, center, ccw);
  return [...arcPieces(a, m, center, ccw, depth + 1), ...arcPieces(m, b, center, ccw, depth + 1)];
}
function halves(p) {
  const m = p.center ? arcPoint(p.a, p.b, p.center, p.ccw) : p.a.add(p.b).scale(HALF);
  return [{ ...p, b: m }, { ...p, a: m }];
}
/** A disc that covers the piece: the chord is a diameter, and an arc of at most a half turn stays inside. */
function cover(p) {
  const d = p.b.sub(p.a);
  return new Disc(p.a.add(p.b).scale(HALF), roundUp(sqrtUp(d.norm().divide(new Rational(4n))).multiply(PAD), 40));
}

/**
 * A path in the plane made of straight pieces and arcs of circles, through
 * exact points.
 *
 *   Contour.polygon(["0", "1", "1 + i"])        closed
 *   Contour.path(["0", "1", "1 + i"])           open
 *   Contour.rectangle(rect([0, 1], [0, 2]))     counter-clockwise
 *   Contour.circle("0", 2)                      counter-clockwise from centre + radius
 *   Contour.arc("0", "1", "i")                  about 0 from 1 to i, counter-clockwise
 */
export class Contour {
  constructor(pieces) {
    if (!pieces.length) throw new RangeError("a contour needs at least one piece");
    for (let k = 1; k < pieces.length; k++) {
      if (!pieces[k].a.equals(pieces[k - 1].b)) throw new RangeError(`the contour is broken: piece ${k} starts at ${pieces[k].a}, the one before ends at ${pieces[k - 1].b}`);
    }
    this.pieces = Object.freeze(pieces.map((p) => Object.freeze({ ...p })));
    Object.freeze(this);
  }
  /** From a Contour, a Rect (its boundary), a Disc (its circle) or an array of points (the closed polygon). */
  static from(x) {
    if (x instanceof Contour) return x;
    if (x instanceof Rect) return Contour.rectangle(x);
    if (x instanceof Disc) return Contour.circle(x.center, x.radius);
    if (Array.isArray(x)) return Contour.polygon(x);
    throw new TypeError("give a Contour, a Rect, a Disc, or an array of points");
  }
  static path(points) {
    const pts = points.map((p) => Complex.from(p));
    if (pts.length < 2) throw new RangeError("a path needs at least two points");
    const pieces = [];
    for (let k = 1; k < pts.length; k++) if (!pts[k].equals(pts[k - 1])) pieces.push(line(pts[k - 1], pts[k]));
    return new Contour(pieces);
  }
  static polygon(points) {
    return Contour.path([...points, points[0]]);
  }
  static rectangle(r) {
    r = Rect.from(r);
    if (r.width.numerator === 0n || r.height.numerator === 0n) throw new RangeError("the rectangle needs a positive width and height");
    return Contour.polygon(r.corners());
  }
  /** The circle, once round. `clockwise: true` reverses it. */
  static circle(center, radius, { clockwise = false } = {}) {
    const c = Complex.from(center);
    const r = Q(radius);
    if (r.numerator <= 0n) throw new RangeError("the circle needs a positive radius");
    const start = c.add(new Complex(r, ZERO));
    return new Contour(arcPieces(start, start, c, !clockwise));
  }
  /** The arc about `center` from one point of a circle to another (the same point twice: the full circle). */
  static arc(center, from, to, { clockwise = false } = {}) {
    const [c, a, b] = [Complex.from(center), Complex.from(from), Complex.from(to)];
    if (!a.sub(c).norm().equals(b.sub(c).norm()) || a.equals(c)) throw new RangeError("the two ends of an arc must be at the same positive distance from its centre");
    return new Contour(arcPieces(a, b, c, !clockwise));
  }
  get start() {
    return this.pieces[0].a;
  }
  get end() {
    return this.pieces[this.pieces.length - 1].b;
  }
  get closed() {
    return this.start.equals(this.end);
  }
  /** This contour followed by another, which must start where this one ends. */
  join(other) {
    return new Contour([...this.pieces, ...Contour.from(other).pieces]);
  }
  /** The same path walked backwards. */
  reversed() {
    return new Contour(this.pieces.map((p) => ({ a: p.b, b: p.a, center: p.center, ccw: !p.ccw })).reverse());
  }
  /** Points along the contour as [x, y] doubles, for drawing: `perPiece` on each piece. */
  sample(perPiece = 8) {
    const out = [this.start.toNumbers()];
    for (const p of this.pieces) {
      for (let k = 1; k <= perPiece; k++) out.push((k === perPiece ? p.b : p.center ? arcPoint(p.a, p.b, p.center, p.ccw, k / perPiece) : p.a.add(p.b.sub(p.a).scale(new Rational(BigInt(k), BigInt(perPiece))))).toNumbers());
    }
    return out;
  }
  /** How many times the closed contour goes round the point w, counter-clockwise positive. */
  winding(w = Complex.ZERO) {
    const c = Complex.from(w);
    return argumentPrinciple((z) => z.sub(c), this, { order: 2 });
  }
  toString() {
    return `${this.closed ? "closed " : ""}contour of ${this.pieces.length} piece${this.pieces.length === 1 ? "" : "s"} from ${this.start}${this.closed ? "" : ` to ${this.end}`}`;
  }
}

// ---- integrals ------------------------------------------------------------------------------

const settings = (opts = {}) => ({ order: opts.order ?? 12, bits: opts.bits ?? 96, maxDepth: opts.maxDepth ?? 30, pieces: { left: opts.maxPieces ?? 20000 } });
/** Each piece tried is paid for: a contour that cannot be settled anywhere would otherwise be halved for ever. */
const spend = (s, p) => {
  if (s.pieces.left-- <= 0) refuse("precision", `the work budget ran out near ${near(p)}: the function is too small or too wild on the contour for this precision (try more bits, or maxPieces)`);
};
const modelOn = (F, D, s) => {
  const z = CTM.z(D, s);
  const out = F(z);
  return out instanceof CTM ? out : z._lift(out);
};
const tol = (t, fallback) => (t === undefined ? fallback : Q(t));
const near = (p) => {
  const [x, y] = p.a.add(p.b).scale(HALF).toNumbers();
  return `${Number(x.toPrecision(6))}${y < 0 ? " - " : " + "}${Number(Math.abs(y).toPrecision(6))}i`;
};

/**
 * The integral of f along a contour, as a rectangle certain to contain it.
 * f must be analytic on the contour (a pole inside a closed contour is fine;
 * one on the contour is refused). Options: tolerance (the width to aim for,
 * default 2^−40), order and bits of the Taylor models, maxDepth.
 */
export function contourIntegral(f, contour, opts = {}) {
  const F = asFunction(f);
  return _integrateMany((z) => [F(z)], contour, opts)[0];
}
/** Several integrands at once over the same contour, sharing the subdivision: G(z) returns an array of models. */
export function _integrateMany(G, contour, opts = {}) {
  const C = Contour.from(contour);
  const s = settings(opts);
  const want = tol(opts.tolerance, new Rational(1n, 1n << 40n)).divide(new Rational(BigInt(C.pieces.length)));
  const lift = (z, m) => (m instanceof CTM ? m : z._lift(m));
  const worst = (rs) => rs.reduce((m, r) => (r.size.compareTo(m) > 0 ? r.size : m), ZERO);
  const piece = (p, depth, goal, before = null) => {
    let got = null;
    let why = null;
    spend(s, p);
    try {
      const z = CTM.z(cover(p), s);
      got = G(z).map((m) => lift(z, m).integral(p.a, p.b));
      if (got.every((r) => r.size.compareTo(goal) <= 0)) return got;
      // half the piece and not half the width: the working precision is the limit, and halving again will not help
      if (before && worst(got).multiply(TWO).compareTo(before) > 0) return got;
    } catch (e) {
      why = refusal(e);
    }
    if (depth >= s.maxDepth) {
      if (got) return got;
      return refuse("singularity", `the function could not be bounded on the contour near ${near(p)}: ${why.message}`);
    }
    const [l, r] = halves(p);
    const g = goal.divide(TWO);
    const was = got ? worst(got) : null;
    const [x, y] = [piece(l, depth + 1, g, was), piece(r, depth + 1, g, was)];
    return x.map((v, k) => v.add(y[k]));
  };
  let total = null;
  for (const p of C.pieces) {
    const v = piece(p, 0, want);
    total = total ? total.map((t, k) => t.add(v[k])) : v;
  }
  return total.map((t) => t.round(s.bits));
}

/**
 * The integral round a closed contour divided by 2πi: by the residue theorem,
 * the sum of the residues of f at the singularities inside (each times the
 * number of times the contour winds round it).
 */
export function residueSum(f, contour, opts = {}) {
  const C = Contour.from(contour);
  if (!C.closed) refuse("bad-argument", "the residue theorem needs a closed contour");
  const s = settings(opts);
  const pi = constants.pi(s.bits);
  const twoPi = I.scale(iv(pi.low, pi.high), TWO);
  // ÷ 2πi: (x + iy)/(2πi) = (y − ix)/2π
  const v = contourIntegral(f, C, opts);
  return new Rect(I.div(v.im, twoPi), I.neg(I.div(v.re, twoPi))).round(s.bits);
}

// ---- the argument principle -------------------------------------------------------------------

/** The change of arg f along one piece, as an interval of radians. */
function argChange(F, p, s, depth = 0) {
  let why = null;
  spend(s, p);
  try {
    const m = modelOn(F, cover(p), s);
    const R = m.range();
    if (!R.containsZero()) {
      const [fa, fb] = [m.at(p.a), m.at(p.b)];
      if (!fa.containsZero() && !fb.containsZero()) {
        const d = arg(fb.div(fa), { bits: 64 });
        // a wide answer is a true one, but many of them add up to more than a turn: ask for better while that is possible
        if (d.value && (depth >= s.maxDepth || I.width(d.value.re).compareTo(SHARP) <= 0)) return d.value.re;
      }
    }
  } catch (e) {
    why = refusal(e);
  }
  if (depth >= s.maxDepth) {
    refuse("zero-on-contour", why ? `the function could not be bounded on the contour near ${near(p)}: ${why.message}` : `the function may vanish on the contour near ${near(p)}`);
  }
  const [l, r] = halves(p);
  return I.add(argChange(F, l, s, depth + 1), argChange(F, r, s, depth + 1));
}

/**
 * The argument principle: the number of times the image f(contour) winds
 * round 0, which for f meromorphic inside a simple closed counter-clockwise
 * contour is (zeros inside) − (poles inside), each counted with its
 * multiplicity. A whole number, proved. Refuses when f may vanish, or cannot
 * be bounded, somewhere on the contour.
 */
export function argumentPrinciple(f, contour, opts = {}) {
  const F = asFunction(f);
  const C = Contour.from(contour);
  if (!C.closed) refuse("bad-argument", "the argument principle needs a closed contour");
  const s = { ...settings(opts), maxDepth: opts.maxDepth ?? 40 };
  let total = iv(ZERO, ZERO);
  for (const p of C.pieces) total = I.add(total, argChange(F, p, s));
  const pi = constants.pi(64);
  const turns = I.div(total, I.scale(iv(pi.low, pi.high), TWO));
  const n = turns.lo.add(turns.hi).divide(TWO).round();
  const k = new Rational(BigInt(n));
  if (!I.has(turns, k) || I.width(turns).compareTo(HALF) > 0) refuse("precision", "the winding number could not be pinned to one whole number");
  return Number(n);
}

// ---- zeros ----------------------------------------------------------------------------------

/**
 * One Newton step on the disc T = (c, ρ): c − f(c)/f′(T). f(z) = f(c) + (z − c)·g
 * with g in the convex hull of f′(T), so every zero in T lies in the result N;
 * and when N is inside T there is exactly one zero in T, a simple one.
 * Returns N, or null when f′ may vanish on T or the function cannot be modelled.
 */
function newton(F, T, s) {
  try {
    const M = modelOn(F, new Disc(T.center, T.radius.multiply(TWO)), s);
    const dF = M.derivative(T.radius).range();
    if (dF.containsZero()) return null;
    return new Disc(T.center).sub(Rect.from(M.at(T.center)).toDisc().div(dF));
  } catch (e) {
    refusal(e);
    return null;
  }
}
const tidy = (D, bits) => {
  const R = Rect.from(D).round(bits);
  return R.toDisc();
};
/** Shrink a proved zero's disc by repeated Newton steps to about 2^−bits. */
function polish(F, N, s, bits) {
  const target = new Rational(1n, 1n << BigInt(bits));
  let cur = N;
  for (let k = 0; k < 20 && cur.radius.compareTo(target) > 0; k++) {
    const t = tidy(cur, bits + 24);
    const T = new Disc(t.center, t.radius.multiply(TWO));
    const next = newton(F, T, s);
    if (!next || !T.contains(next) || next.radius.multiply(TWO).compareTo(cur.radius) > 0) break;
    cur = next;
  }
  return cur;
}

/**
 * The zeros of an analytic f in a rectangle, isolated and counted.
 *
 *   zeros("exp(z) - 3z", rect([-1, 3], [-2, 2]))
 *   → { zeros: [{ enclosure, multiplicity, simple, isolating, boundary }], count, complete, unresolved }
 *
 * Each entry's `enclosure` (a Rect) holds zeros of total multiplicity
 * `multiplicity`. When `simple` is true it is one simple zero, the only zero
 * in the disc `isolating`. Otherwise the zeros in the box could not be
 * separated at the tolerance: a multiple zero, or several closer together
 * than `tolerance`, and `multiplicity` is their proven total.
 *
 * `complete` is true when every part of the rectangle was settled, and then
 * `count` is the proven number of zeros in it (entries flagged `boundary`
 * straddle the edge and are not included in `count`). `unresolved` lists
 * the parts that could not be settled, with reasons: a singularity, or the
 * work budget.
 *
 * Options: tolerance (how small a box an unseparated cluster is narrowed to,
 * default 2^−20), digits of each simple zero via bits (default 64), order and
 * modelBits of the Taylor models, maxCells (the work budget, default 20000).
 */
export function zeros(f, region, opts = {}) {
  const F = asFunction(f);
  const R = Rect.from(region);
  if (R.width.numerator === 0n || R.height.numerator === 0n) throw new RangeError("the region needs a positive width and height");
  const bits = opts.bits ?? 64;
  const s = { order: opts.order ?? 10, bits: Math.max(opts.modelBits ?? 96, bits + 40), maxDepth: 40 };
  const tolerance = tol(opts.tolerance, new Rational(1n, 1n << 20n));
  const maxCells = opts.maxCells ?? 20000;
  const FIVE_FOURTHS = new Rational(5n, 4n);
  // near a zero of order k the function is of size radius^k: fixed-point models of small discs need bits to match
  const deep = (radius, k) => {
    const lg = radius.denominator.toString(2).length - radius.numerator.toString(2).length;
    return { ...s, bits: Math.min(4000, s.bits + k * Math.max(0, lg)) };
  };

  // near-square starting cells
  const [w, h] = [R.width, R.height];
  const ratio = w.compareTo(h) >= 0 ? w.divide(h) : h.divide(w);
  const n = Math.max(1, Math.min(64, Math.round(ratio.toNumber())));
  let cells = [];
  for (let k = 0; k < n; k++) {
    const part = (side) => iv(side.lo.add(I.width(side).multiply(new Rational(BigInt(k), BigInt(n)))), side.lo.add(I.width(side).multiply(new Rational(BigInt(k + 1), BigInt(n)))));
    cells.push(w.compareTo(h) >= 0 ? new Rect(part(R.re), R.im) : new Rect(R.re, part(R.im)));
  }

  const found = []; // { iso, enc, k }: the k zeros in the disc iso are all in the disc enc
  const leftover = [];
  const unresolved = [];
  let used = 0;
  const keep = (iso, enc, k) => {
    let doubt = false;
    for (const e of found) {
      if (e.iso.contains(enc) || iso.contains(e.enc)) return; // the same zeros, found from a neighbouring cell
      if (e.enc.intersects(enc)) doubt = true;
    }
    found.push({ iso, enc, k, doubt });
  };
  const record = (iso, N) => keep(iso, polish(F, N, s, bits), 1);
  /**
   * A disc proved to hold k ≥ 2 zeros: follow them inward, re-centring on
   * their centroid, for as long as a smaller disc still holds all k.
   */
  const narrow = (m, k) => {
    let cur = m;
    // past the tolerance too, while it still works: k zeros can be told apart from one zero of order k only
    // down to about the k-th root of the working precision
    const floor = new Rational(1n, 1n << BigInt(Math.max(20, Math.floor(bits / k))));
    for (let step = 0; step < 200 && cur.radius.compareTo(floor) > 0; step++) {
      // the centroid of the zeros, from the two top terms: −a_(k−1)/(k·a_k), in units of the radius
      const shift = cur._c(cur.a[k - 1] ?? [0n, 0n]).div(cur._c(cur.a[k])).scale(new Rational(-1n, BigInt(k)));
      const c = Rect.from(cur.center.add(shift.scale(cur.radius))).round(cur.bits).center;
      let next = null;
      for (const f of [new Rational(1n, 4n), HALF]) {
        const D = new Disc(c, roundUp(cur.radius.multiply(f), 40));
        if (!cur.disc.contains(D)) continue;
        try {
          const t = modelOn(F, D, deep(D.radius, k + 1));
          if (t.zeroCount() === k) next = t;
        } catch (e) {
          refusal(e);
        }
        if (next) break;
      }
      if (!next) break;
      cur = next;
    }
    return cur.disc;
  };
  while (cells.length) {
    const next = [];
    for (let i = 0; i < cells.length; i++) {
      const cell = cells[i];
      if (used >= maxCells) {
        unresolved.push({ region: Rect.hull(...cells.slice(i), ...next), reason: `the work budget of ${maxCells} cells ran out; the function may be identically zero here, or the region may need to be smaller` });
        next.length = 0;
        break;
      }
      const D = cell.toDisc();
      if (found.some((e) => e.iso.contains(D))) continue; // inside a disc already known to hold exactly one zero
      used++;
      const small = cell.size.compareTo(tolerance) <= 0;
      let m = null;
      let why = null;
      try {
        m = modelOn(F, new Disc(D.center, roundUp(D.radius, 40)), deep(D.radius, 3));
      } catch (e) {
        why = refusal(e);
      }
      if (m && !m.range().containsZero()) continue; // no zero here
      const k = m ? m.zeroCount() : null;
      if (k === 0) continue;
      if (k >= 2) {
        const inner = narrow(m, k);
        if (inner.radius.compareTo(tolerance) <= 0) {
          keep(m.disc, inner, k);
          continue;
        }
        if (inner.radius.multiply(TWO).compareTo(m.radius) <= 0) {
          // every zero of this cell is in the smaller disc: carry on from there
          next.push(Rect.from(inner));
          continue;
        }
      }
      if (m) {
        const T = new Disc(D.center, roundUp(D.radius.multiply(FIVE_FOURTHS), 40));
        const N = newton(F, T, s);
        if (N && T.contains(N)) {
          record(T, N);
          continue;
        }
      }
      if (small) {
        if (m) leftover.push(cell);
        else unresolved.push({ region: cell, reason: `the function could not be bounded here: ${why.message}` });
        continue;
      }
      next.push(...cell.split());
    }
    cells = next;
  }

  // what is left are small cells where a zero could be neither excluded nor shown simple: group and count
  const clusters = [];
  if (leftover.length > 400) {
    unresolved.push({ region: Rect.hull(...leftover), reason: "too many cells could not be settled; the function may be identically zero" });
  } else {
    const group = leftover.map((_, i) => i);
    const top = (i) => (group[i] === i ? i : (group[i] = top(group[i])));
    const fat = leftover.map((c) => grow(c, c.size.divide(TWO)));
    for (let i = 0; i < leftover.length; i++) for (let j = 0; j < i; j++) if (fat[i].intersects(fat[j])) group[top(i)] = top(j);
    const boxes = new Map();
    leftover.forEach((c, i) => boxes.set(top(i), boxes.has(top(i)) ? boxes.get(top(i)).hull(c) : c));
    for (const box of boxes.values()) {
      const B = grow(box, tolerance);
      try {
        const sB = deep(B.size, 6);
        modelOn(F, B.toDisc(), sB); // analytic on the whole box
        const k = argumentPrinciple(F, Contour.rectangle(B), { order: sB.order, bits: sB.bits, maxPieces: 4000 });
        if (k > 0) {
          // k zeros in a small box: centre a disc on their centroid and follow them inward
          let box = B;
          try {
            const D0 = B.toDisc();
            const m0 = modelOn(F, new Disc(D0.center, D0.radius.multiply(TWO)), sB);
            const shift = m0._c(m0.a[k - 1] ?? [0n, 0n]).div(m0._c(m0.a[k])).scale(new Rational(-1n, BigInt(k)));
            const c1 = Rect.from(m0.center.add(shift.scale(m0.radius))).round(sB.bits).center;
            const D1 = new Disc(c1, roundUp(sqrtUp(c1.sub(D0.center).norm()).add(D0.radius), 40));
            const m1 = modelOn(F, D1, deep(D1.radius, k + 1));
            if (m1.zeroCount() === k) box = Rect.from(narrow(m1, k)).round(sB.bits);
          } catch (e) {
            refusal(e);
          }
          clusters.push({ box, k });
        }
        else if (k < 0) unresolved.push({ region: B, reason: "the count came out negative: the function is not analytic here" });
      } catch (e) {
        unresolved.push({ region: B, reason: refusal(e).message });
      }
    }
  }

  const out = [];
  for (const e of found) {
    const enclosure = Rect.from(e.enc).round(bits + 8);
    if (!R.intersects(enclosure)) continue;
    out.push({ enclosure, multiplicity: e.k, simple: e.k === 1, isolating: e.iso, boundary: !R.contains(enclosure), ...(e.doubt ? { possibleDuplicate: true } : {}) });
  }
  for (const c of clusters) {
    if (!R.intersects(c.box)) continue;
    // a zero already listed that lies in the box is part of the box's count
    const insideBox = out.filter((z) => c.box.contains(z.enclosure)).reduce((t, z) => t + z.multiplicity, 0);
    const straddling = out.some((z) => c.box.intersects(z.enclosure) && !c.box.contains(z.enclosure));
    if (straddling) {
      unresolved.push({ region: c.box, reason: "a cluster of zeros overlaps a zero already found" });
      continue;
    }
    if (c.k - insideBox > 0) out.push({ enclosure: c.box, multiplicity: c.k - insideBox, simple: false, isolating: null, boundary: !R.contains(c.box) });
  }
  const key = (z) => z.enclosure.center;
  out.sort((a, b) => key(a).re.compareTo(key(b).re) || key(a).im.compareTo(key(b).im));
  const complete = unresolved.length === 0;
  const count = out.filter((z) => !z.boundary).reduce((t, z) => t + z.multiplicity, 0);
  return { zeros: out, count, complete, unresolved };
}
const grow = (r, by) => new Rect(iv(r.re.lo.subtract(by), r.re.hi.add(by)), iv(r.im.lo.subtract(by), r.im.hi.add(by)));

/**
 * The number of zeros of an analytic f inside a closed contour (or a Rect or
 * Disc), counted with multiplicity. The same number as argumentPrinciple;
 * the name says what it means when f has no poles inside.
 */
export const countZeros = (f, contour, opts) => argumentPrinciple(f, contour, opts);
