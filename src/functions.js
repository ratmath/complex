// Polar form and the elementary functions, on enclosures.
//
// Single-valued functions (exp, sin, cos, sinh, cosh, tan, tanh) take a
// Complex, a Rect or a Disc and return a Rect certain to contain every value.
//
// Multi-valued functions (log, sqrt, roots, powers, the inverse trigonometric
// and hyperbolic functions, arg) return a Branched result, because the
// principal branch is a choice and the choice must not be made silently: when
// the input touches a branch cut the result says so and carries the value on
// each side. A branch can be named ({ branch: k }), a side chosen
// ({ side: "above" }), or a value carried continuously along a path
// (continueAlong).

import { Rational, RationalInterval } from "@ratmath/core";
import { enclose, constants } from "@ratmath/analytic";
import { Refinement } from "@ratmath/core/refine";
import { Complex, Q, ZERO, ONE, TWO, refuse } from "./complex.js";
import { Rect, Disc, I, iv, sqrtUp, sqrtDown } from "./enclosure.js";

const HALF = new Rational(1n, 2n);
const NEG1 = new Rational(-1n);
const EXTRA = 16;

const F = (name, a, bits) => {
  const r = enclose[name](I.isPoint(a) ? a.lo : I.toCore(a), bits);
  return iv(r.low, r.high);
};
const piIv = (bits) => {
  const p = constants.pi(bits);
  return iv(p.low, p.high);
};
const pt = (q) => iv(q, q);
/**
 * atan at an exact point. A point whose denominator is longer than the working
 * precision is first widened to a dyadic interval of that precision: atan′ ≤ 1,
 * so nothing is lost, and analytic's series is never handed an oversized fraction.
 */
function atanAt(t, bits) {
  if (t.denominator.toString(2).length <= bits) return F("atan", pt(t), bits);
  const scale = 1n << BigInt(bits);
  const n = t.numerator * scale;
  let f = n / t.denominator;
  if (n % t.denominator !== 0n && n < 0n) f -= 1n;
  return F("atan", iv(new Rational(f, scale), new Rational(f + 1n, scale)), bits);
}
const toRect = (z) => Rect.from(z);
const exactInput = (z) => (z instanceof Rect || z instanceof Disc ? z.isPoint : true);
/** Working precision and the final rounding: exact inputs are computed with guard bits, then rounded outward. */
const plan = (z, opts) => {
  const bits = typeof opts === "number" ? opts : opts?.bits ?? 64;
  return { bits, work: bits + EXTRA, finish: (r) => (exactInput(z) ? r.round(bits + 4) : r.round(bits + EXTRA)) };
};

// ---- modulus and argument ---------------------------------------------------------------

/** |z| as a Core RationalInterval; a single point when the modulus is rational. */
export function modulus(z, bits = 64) {
  const r = toRect(z);
  const n = r.normRange();
  return new RationalInterval(sqrtDown(n.lo, bits), sqrtUp(n.hi, bits));
}

/** atan2 at an exact point; `onAxis` says which value the negative real axis gets: π ("above") or −π ("below"). */
function argPoint(x, y, bits, onAxis = "above") {
  const sx = x.numerator;
  const sy = y.numerator;
  if (sx === 0n && sy === 0n) refuse("branch-point", "the argument of 0 is undefined");
  const pi = piIv(bits);
  const half = I.scale(pi, HALF);
  if (sx === 0n) return sy > 0n ? half : I.neg(half);
  if (y.abs().compareTo(x.abs()) > 0) {
    // nearer the imaginary axis: arg = ±π/2 − atan(x/y), which keeps the arctangent's argument small
    const a = atanAt(x.divide(y), bits);
    return sy > 0n ? I.sub(half, a) : I.sub(I.neg(half), a);
  }
  const a = atanAt(y.divide(x), bits);
  if (sx > 0n) return a;
  return sy > 0n || (sy === 0n && onAxis === "above") ? I.add(a, pi) : I.sub(a, pi);
}

const touchesNegativeAxis = (r) => r.re.lo.numerator < 0n && I.hasZero(r.im);
const clip = (r, side) => {
  switch (side) {
    case "above": return new Rect(r.re, iv(r.im.lo.numerator < 0n ? ZERO : r.im.lo, r.im.hi));
    case "below": return new Rect(r.re, iv(r.im.lo, r.im.hi.numerator > 0n ? ZERO : r.im.hi));
    case "right": return new Rect(iv(r.re.lo.numerator < 0n ? ZERO : r.re.lo, r.re.hi), r.im);
    case "left": return new Rect(iv(r.re.lo, r.re.hi.numerator > 0n ? ZERO : r.re.hi), r.im);
  }
  return r;
};

/**
 * The range of arg over a rectangle that does not contain 0, as { lo, hi }.
 * With side "above" the rectangle is taken to lie in Im ≥ 0 and the negative
 * real axis has argument π; with "below", Im ≤ 0 and −π. Without a side the
 * rectangle must not touch the negative real axis.
 */
function argRange(r, bits, side = null) {
  if (r.containsZero()) refuse("branch-point", `the enclosure ${r} contains 0, where the argument is undefined`);
  if (side === null && touchesNegativeAxis(r)) refuse("branch-cut", `the enclosure ${r} touches the negative real axis`);
  // arg has no interior extremum and is monotone along each edge, so the corners bound it
  return r.corners().map((c) => argPoint(c.re, c.im, bits, side ?? "above")).reduce(I.hull);
}

// ---- results with branches ----------------------------------------------------------------

/**
 * The result of a multi-valued function.
 *
 *   value       the enclosure (a Rect), or null when the input straddles the cut
 *   straddles   true when the input touches the branch cut
 *   cut         the cut of this function, in words
 *   candidates  when it straddles: [{ side, value }] for each side of the cut
 *   convention  when an exact input lies on the cut: the side the usual convention takes, and its value
 *   branch      which branch `value` is on: "principal", or the number given
 */
export class Branched {
  constructor(fields) {
    Object.assign(this, { value: null, straddles: false, candidates: [], convention: null, branch: "principal" }, fields);
    Object.freeze(this.candidates);
    Object.freeze(this);
  }
  /** The enclosure, or an error that says why there is not one. */
  get rect() {
    if (this.value) return this.value;
    return refuse("branch-cut", `${this.name}: the input touches the branch cut (${this.cut}); choose { side: ${this.candidates.map((c) => `"${c.side}"`).join(" or ")} }`);
  }
  /** A rectangle around every candidate (the function's whole jump across the cut is inside it). */
  get hull() {
    return this.value ?? this.candidates.map((c) => c.value).reduce((a, b) => a.hull(b));
  }
  toString(digits) {
    if (this.value) return `${this.value.toString(digits)}${this.branch === "principal" ? "" : `  (branch ${this.branch})`}`;
    const parts = this.candidates.map((c) => `${c.side}: ${c.value.toString(digits)}`).join(";  ");
    return `on the cut (${this.cut}).  ${parts}${this.convention ? `.  Convention: ${this.convention.side}` : ""}`;
  }
}

/**
 * Shared shape of every function with a cut on an axis.
 *   touches(r)        does the rectangle meet the cut?
 *   sides             the two sides' names
 *   evalSide(r, side) the value on the rectangle clipped to that side (side null: away from the cut)
 *   convention(z)     for an exact point on the cut, the side the usual convention takes
 */
function withCut(name, z, opts, { cut, touches, sides, evalSide, convention, shift, infer }) {
  const p = plan(z, opts);
  const o = typeof opts === "object" && opts ? opts : {};
  const r = toRect(z);
  const done = (v) => p.finish(shift ? shift(v, o.branch ?? 0, p.work) : v);
  const branch = shift && o.branch ? o.branch : "principal";
  if (o.branch && !shift) refuse("bad-option", `${name} has no numbered branches; use { side } on its cut`);
  if (!touches(r)) {
    // away from the cut, but still on one side of the axis: the inner functions are told which
    let side = null;
    if (infer === "real") side = r.im.lo.numerator >= 0n ? "above" : r.im.hi.numerator <= 0n ? "below" : null;
    if (infer === "imag") side = r.re.lo.numerator >= 0n ? "right" : r.re.hi.numerator <= 0n ? "left" : null;
    return new Branched({ name, cut, branch, value: done(evalSide(r, side, p.work)) });
  }
  if (o.side) {
    if (!sides.includes(o.side)) refuse("bad-option", `${name}: side must be ${sides.map((s) => `"${s}"`).join(" or ")}`);
    return new Branched({ name, cut, branch, value: done(evalSide(clip(r, o.side), o.side, p.work)) });
  }
  const candidates = sides.map((side) => ({ side, value: done(evalSide(clip(r, side), side, p.work)) }));
  const conv = r.isPoint && convention ? candidates.find((c) => c.side === convention(r.center)) : null;
  return new Branched({ name, cut, branch, straddles: true, candidates, convention: conv ?? null });
}

const REAL_SIDES = ["above", "below"];
const IMAG_SIDES = ["right", "left"];
const NEG_AXIS = "the negative real axis";

// ---- exp and friends (single-valued) ----------------------------------------------------------

function expRect(r, bits) {
  const e = F("exp", r.re, bits);
  return new Rect(I.mul(e, F("cos", r.im, bits)), I.mul(e, F("sin", r.im, bits)));
}
function sinRect(r, bits) {
  // sin(x + iy) = sin x cosh y + i cos x sinh y
  return new Rect(I.mul(F("sin", r.re, bits), F("cosh", r.im, bits)), I.mul(F("cos", r.re, bits), F("sinh", r.im, bits)));
}
function cosRect(r, bits) {
  // cos(x + iy) = cos x cosh y − i sin x sinh y
  return new Rect(I.mul(F("cos", r.re, bits), F("cosh", r.im, bits)), I.neg(I.mul(F("sin", r.re, bits), F("sinh", r.im, bits))));
}
// Taylor models (taylor.js) register here, so that exp(model) and the rest work on them.
let modelHook = null;
export const _setModelHook = (h) => (modelHook = h);
const single = (f, name) => (z, opts) => {
  if (modelHook?.is(z)) return modelHook.apply(name, [z]);
  const p = plan(z, opts);
  return p.finish(f(toRect(z), p.work));
};
/** e^z. */
export const exp = single(expRect, "exp");
export const sin = single(sinRect, "sin");
export const cos = single(cosRect, "cos");
/** sinh z = −i sin(iz). */
export const sinh = single((r, b) => sinRect(r.timesI(), b).timesI().neg(), "sinh");
/** cosh z = cos(iz). */
export const cosh = single((r, b) => cosRect(r.timesI(), b), "cosh");
const quotient = (num, den, name) => (r, b) => {
  const d = den(r, b);
  if (d.containsZero()) refuse("pole", `${name}: the enclosure ${r} may contain a pole`);
  return num(r, b).div(d);
};
export const tan = single(quotient(sinRect, cosRect, "tan"), "tan");
export const tanh = single((r, b) => quotient(sinRect, cosRect, "tanh")(r.timesI(), b).timesI().neg(), "tanh");

// ---- log, roots, powers ---------------------------------------------------------------------

function logRect(r, bits, side) {
  if (r.containsZero()) refuse("branch-point", `log: the enclosure ${r} contains 0`);
  const n = r.normRange();
  return new Rect(I.scale(F("log", n, bits), HALF), argRange(r, bits, side));
}
const addTurns = (v, k, bits) => (k ? v.add(new Rect(pt(ZERO), I.scale(piIv(bits), new Rational(BigInt(2 * k))))) : v);

/**
 * The logarithm. The principal branch has −π < Im ≤ π with its cut along the
 * negative real axis; { branch: k } adds 2πik.
 */
export function log(z, opts) {
  if (modelHook?.is(z)) return modelHook.apply("log", [z, opts]);
  return withCut("log", z, opts, { cut: NEG_AXIS, touches: touchesNegativeAxis, sides: REAL_SIDES, evalSide: (r, side, bits) => logRect(r, bits, side), convention: () => "above", shift: addTurns });
}

function rootRect(r, n, bits, side) {
  if (r.isPoint && r.center.isZero()) return r;
  if (r.containsZero()) refuse("branch-point", `a root of the enclosure ${r}, which contains its branch point 0`);
  const nr = r.normRange();
  // |z|^(1/n) = exp(log|z|² / 2n); for square roots use the square root twice
  const m = n === 2 ? F("sqrt", F("sqrt", nr, bits), bits) : F("exp", I.scale(F("log", nr, bits), new Rational(1n, BigInt(2 * n))), bits);
  const th = I.scale(argRange(r, bits, side), new Rational(1n, BigInt(n)));
  return new Rect(I.mul(m, F("cos", th, bits)), I.mul(m, F("sin", th, bits)));
}
const turnRoot = (n) => (v, k, bits) => {
  if (!k) return v;
  const th = I.scale(piIv(bits), new Rational(BigInt(2 * k), BigInt(n)));
  return v.mul(new Rect(F("cos", th, bits), F("sin", th, bits)));
};

/**
 * The n-th root. The principal root has −π/n < arg ≤ π/n; { branch: k }
 * multiplies it by e^(2πik/n). Exact when the root is a complex rational.
 */
export function root(z, n, opts) {
  if (modelHook?.is(z)) return modelHook.apply("root", [z, n, opts]);
  if (!Number.isInteger(n) || n < 1) throw new RangeError("root needs a positive whole number n");
  const o = typeof opts === "object" && opts ? opts : {};
  if (n === 2 && exactInput(z) && o.side) {
    const c = toRect(z).center;
    const w = c.exactSqrt();
    if (w && c.im.equals(ZERO) && c.re.numerator < 0n) {
      const up = w[0].im.numerator > 0n ? w[0] : w[1];
      return new Branched({ name: "sqrt", cut: NEG_AXIS, branch: o.branch ? o.branch : "principal", value: Rect.from(((o.side === "above") !== Boolean(o.branch % 2)) ? up : up.neg()) });
    }
  }
  if (n === 2 && exactInput(z) && !o.side) {
    // an exact square root, when there is one
    const c = toRect(z).center;
    const w = c.exactSqrt();
    const onCut = c.im.equals(ZERO) && c.re.numerator < 0n;
    const branch = o.branch ? o.branch : "principal";
    const pick = (v) => Rect.from(o.branch % 2 ? v.neg() : v);
    if (w && !onCut) return new Branched({ name: "sqrt", cut: NEG_AXIS, branch, value: pick(w[0]) });
    if (w) {
      // a negative real number with a rational root: ±(a positive multiple of i), one on each side of the cut
      const up = w[0].im.numerator > 0n ? w[0] : w[1];
      const candidates = [{ side: "above", value: pick(up) }, { side: "below", value: pick(up.neg()) }];
      return new Branched({ name: "sqrt", cut: NEG_AXIS, branch, straddles: true, candidates, convention: candidates[0] });
    }
  }
  return withCut(n === 2 ? "sqrt" : `root(${n})`, z, opts, {
    cut: NEG_AXIS, touches: (r) => n > 1 && touchesNegativeAxis(r), sides: REAL_SIDES,
    evalSide: (r, side, bits) => rootRect(r, n, bits, side), convention: () => "above", shift: turnRoot(n),
  });
}
/** The square root: root(z, 2). sqrt("3 + 4i") is exactly 2 + i. */
export const sqrt = (z, opts) => root(z, 2, opts);

/** All n values of the n-th root of an exact z, counter-clockwise from the principal one. */
export function roots(z, n, opts) {
  const bits = typeof opts === "number" ? opts : opts?.bits ?? 64;
  const c = Complex.from(z);
  if (c.isZero()) return [Rect.from(Complex.ZERO)];
  return Array.from({ length: n }, (_, k) => root(c, n, { bits, branch: k, side: "above" }).rect);
}

/** The n-th roots of unity, e^(2πik/n) for k = 0 … n − 1; the ones that are complex rationals are exact points. */
export function rootsOfUnity(n, bits = 64) {
  const exact = { 0: Complex.ONE, 1: Complex.I, 2: Complex.ONE.neg(), 3: Complex.I.neg() };
  return Array.from({ length: n }, (_, k) => {
    if ((4 * k) % n === 0) return Rect.from(exact[(4 * k) / n]);
    const th = I.scale(piIv(bits + EXTRA), new Rational(BigInt(2 * k), BigInt(n)));
    return new Rect(F("cos", th, bits + EXTRA), F("sin", th, bits + EXTRA)).round(bits + 4);
  });
}

/**
 * z^w = exp(w·log z). A whole-number power of an exact z is exact. Otherwise
 * the branch of the logarithm matters, and the result reports it the way log does.
 */
export function pow(z, w, opts) {
  if (modelHook?.is(z) || modelHook?.is(w)) return modelHook.apply("pow", [z, w, opts]);
  const o = typeof opts === "object" && opts ? opts : {};
  const wr = toRect(w);
  if (wr.isPoint && wr.center.isReal() && wr.center.re.denominator === 1n) {
    const k = Number(wr.center.re.numerator);
    const r = toRect(z);
    const p = plan(z, opts);
    return new Branched({ name: "pow", cut: "none: a whole-number power is single-valued", value: r.isPoint ? Rect.from(r.center.pow(k)) : p.finish(r.pow(k)) });
  }
  if (wr.isPoint && wr.center.isReal() && wr.center.re.denominator <= 64n && wr.center.re.numerator >= -64n && wr.center.re.numerator <= 64n) {
    // z^(p/q) is the principal q-th root to the power p
    const e = wr.center.re;
    const r = root(z, Number(e.denominator), opts);
    const up = (v) => (v.isPoint ? Rect.from(v.center.pow(Number(e.numerator))) : v.pow(Number(e.numerator)).round((typeof opts === "number" ? opts : opts?.bits ?? 64) + 4));
    return new Branched({ name: "pow", cut: NEG_AXIS, branch: r.branch, straddles: r.straddles, value: r.value && up(r.value), candidates: r.candidates.map((c) => ({ side: c.side, value: up(c.value) })), convention: r.convention && { side: r.convention.side, value: up(r.convention.value) } });
  }
  const p = plan(z, opts);
  const l = log(z, { bits: p.work, branch: o.branch, side: o.side });
  const f = (v) => p.finish(expRect(wr.mul(v), p.work));
  if (!l.straddles) return new Branched({ name: "pow", cut: NEG_AXIS, branch: l.branch, value: f(l.value) });
  const candidates = l.candidates.map((c) => ({ side: c.side, value: f(c.value) }));
  return new Branched({ name: "pow", cut: NEG_AXIS, branch: l.branch, straddles: true, candidates, convention: l.convention ? candidates.find((c) => c.side === l.convention.side) : null });
}

/**
 * The argument. Principal values lie in (−π, π]; the result is a Branched
 * whose rectangles have real part arg and imaginary part 0. Use
 * argument() for the plain interval.
 */
export function arg(z, opts) {
  const lift = (r, side, bits) => new Rect(argRange(r, bits, side), pt(ZERO));
  return withCut("arg", z, opts, { cut: NEG_AXIS, touches: touchesNegativeAxis, sides: REAL_SIDES, evalSide: lift, convention: () => "above", shift: (v, k, bits) => (k ? v.add(new Rect(I.scale(piIv(bits), new Rational(BigInt(2 * k))), pt(ZERO))) : v) });
}

/** The principal argument of an exact z as a Core RationalInterval in (−π, π]. */
export function argument(z, bits = 64) {
  const c = Complex.from(z);
  return I.toCore(I.round(argPoint(c.re, c.im, bits + EXTRA), bits + 4));
}

/**
 * The exact fraction of a full turn when it is rational, else null. For a
 * complex rational that happens only on the axes and the diagonals: the turn
 * is a multiple of 1/8.
 */
export function exactTurn(z) {
  const c = Complex.from(z);
  if (c.isZero()) return null;
  const sx = c.re.numerator > 0n ? 1 : c.re.numerator < 0n ? -1 : 0;
  const sy = c.im.numerator > 0n ? 1 : c.im.numerator < 0n ? -1 : 0;
  if (sx !== 0 && sy !== 0 && !c.re.abs().equals(c.im.abs())) return null;
  const eighths = { "1,0": 0, "1,1": 1, "0,1": 2, "-1,1": 3, "-1,0": 4, "-1,-1": -3, "0,-1": -2, "1,-1": -1 }[`${sx},${sy}`];
  return new Rational(BigInt(eighths), 8n);
}

/**
 * Polar form of an exact z: { modulus, argument } as Core RationalIntervals,
 * with what is known exactly: normSquared, exactModulus (or null) and
 * exactTurn (the argument as a fraction of a full turn, or null).
 */
export function polar(z, bits = 64) {
  const c = Complex.from(z);
  return { modulus: modulus(c, bits), argument: argument(c, bits), normSquared: c.norm(), exactModulus: c.exactModulus(), exactTurn: exactTurn(c) };
}

/**
 * r·e^(iθ) as a Rect. θ is in radians (a Rational or a RationalInterval), or
 * { turns: q } for the exact angle 2πq, which gives exact points on the axes.
 */
export function fromPolar(r, theta, bits = 64) {
  const m = I.of(r);
  if (theta && typeof theta === "object" && "turns" in theta) {
    const t = Q(theta.turns);
    const quarter = t.multiply(new Rational(4n));
    if (quarter.denominator === 1n) {
      const k = Number(((quarter.numerator % 4n) + 4n) % 4n);
      return new Rect(pt(ZERO), pt(ZERO)).add(Rect.from([Complex.ONE, Complex.I, Complex.ONE.neg(), Complex.I.neg()][k]).mul(new Rect(m, pt(ZERO))));
    }
    const th = I.scale(piIv(bits + EXTRA), t.multiply(TWO));
    return new Rect(I.mul(m, F("cos", th, bits + EXTRA)), I.mul(m, F("sin", th, bits + EXTRA))).round(bits + 4);
  }
  const th = I.of(theta);
  return new Rect(I.mul(m, F("cos", th, bits + EXTRA)), I.mul(m, F("sin", th, bits + EXTRA))).round(bits + 4);
}

// ---- inverse functions ------------------------------------------------------------------------

const one = Rect.from(Complex.ONE);
const flip = (side) => (side === "above" ? "below" : side === "below" ? "above" : null);
// An inner value known to lie on one side of the real axis: an enclosure of it that pokes across is clipped back.
const sqrtSide = (r, bits, side) => (touchesNegativeAxis(r) ? rootRect(clip(r, side ?? refuseCut(r)), 2, bits, side) : rootRect(r, 2, bits, null));
const logSide = (r, bits, side) => (touchesNegativeAxis(r) ? logRect(clip(r, side ?? refuseCut(r)), bits, side) : logRect(r, bits, null));
const refuseCut = (r) => refuse("precision", `an intermediate enclosure ${r} is too wide to place relative to a branch cut; use a smaller input or more bits`);

/** asin on a rectangle on one side of the real axis (side null: away from the cut). */
function asinRect(r, bits, side) {
  // asin z = −i log(iz + S), S = sqrt(1 − z)·sqrt(1 + z). Above the axis 1 − z is below it and 1 + z above.
  const S = sqrtSide(one.sub(r), bits, flip(side)).mul(sqrtSide(one.add(r), bits, side));
  const iz = r.timesI();
  // (iz + S)(S − iz) = 1: use whichever sum does not cancel
  const l = r.im.lo.numerator >= 0n ? logSide(S.sub(iz), bits, null).neg() : logSide(iz.add(S), bits, null);
  return l.timesI().neg();
}
const realCut = (r) => I.hasZero(r.im) && (r.re.lo.compareTo(NEG1) < 0 || r.re.hi.compareTo(ONE) > 0);
const outsideUnit = "the real axis outside [-1, 1]";
// On a cut, "convention" is what the defining formula gives with principal log and sqrt. Other systems choose differently there.
const formulaReal = (c) => (c.re.numerator > 0n ? "below" : "above");

/** The inverse sine. Its cuts are the real axis left of −1 and right of 1. */
export function asin(z, opts) {
  return withCut("asin", z, opts, { cut: outsideUnit, touches: realCut, sides: REAL_SIDES, evalSide: (r, side, bits) => asinRect(r, bits, side), convention: formulaReal, infer: "real" });
}
/** The inverse cosine, π/2 − asin z. */
export function acos(z, opts) {
  const half = (bits) => new Rect(I.scale(piIv(bits), HALF), pt(ZERO));
  return withCut("acos", z, opts, { cut: outsideUnit, touches: realCut, sides: REAL_SIDES, evalSide: (r, side, bits) => half(bits).sub(asinRect(r, bits, side)), convention: formulaReal, infer: "real" });
}
function atanhRect(r, bits, side) {
  // atanh z = (log(1 + z) − log(1 − z)) / 2
  return logSide(one.add(r), bits, side).sub(logSide(one.sub(r), bits, flip(side))).scale(HALF);
}
/** The inverse hyperbolic tangent. Its cuts are the real axis left of −1 and right of 1; ±1 are singular. */
export function atanh(z, opts) {
  return withCut("atanh", z, opts, { cut: outsideUnit, touches: realCut, sides: REAL_SIDES, evalSide: (r, side, bits) => atanhRect(r, bits, side), convention: formulaReal, infer: "real" });
}
const imagCut = (r) => I.hasZero(r.re) && (r.im.lo.compareTo(NEG1) < 0 || r.im.hi.compareTo(ONE) > 0);
const outsideUnitImag = "the imaginary axis outside [-i, i]";
// multiplying by i turns "right of the imaginary axis" into "above the real axis"
const turned = (side) => (side === "right" ? "above" : side === "left" ? "below" : null);
const formulaImag = (c) => (c.im.numerator > 0n ? "right" : "left");
/** The inverse tangent, −i·atanh(iz). Its cuts are the imaginary axis above i and below −i. */
export function atan(z, opts) {
  return withCut("atan", z, opts, { cut: outsideUnitImag, touches: imagCut, sides: IMAG_SIDES, evalSide: (r, side, bits) => atanhRect(r.timesI(), bits, turned(side)).timesI().neg(), convention: formulaImag, infer: "imag" });
}
/** The inverse hyperbolic sine, −i·asin(iz). Its cuts are the imaginary axis above i and below −i. */
export function asinh(z, opts) {
  return withCut("asinh", z, opts, { cut: outsideUnitImag, touches: imagCut, sides: IMAG_SIDES, evalSide: (r, side, bits) => asinRect(r.timesI(), bits, turned(side)).timesI().neg(), convention: formulaImag, infer: "imag" });
}
/** The inverse hyperbolic cosine, log(z + sqrt(z − 1)·sqrt(z + 1)). Its cut is the real axis left of 1. */
export function acosh(z, opts) {
  const f = (r, side, bits) => logSide(r.add(sqrtSide(r.sub(one), bits, side).mul(sqrtSide(r.add(one), bits, side))), bits, side);
  return withCut("acosh", z, opts, { cut: "the real axis left of 1", touches: (r) => I.hasZero(r.im) && r.re.lo.compareTo(ONE) < 0, sides: REAL_SIDES, evalSide: f, convention: () => "above", infer: "real" });
}

// ---- continuation along a path ------------------------------------------------------------------

/**
 * How many times the closed polygon through the points winds around w
 * (default 0), counter-clockwise positive. Exact. The polygon must not pass through w.
 */
export function winding(path, w = Complex.ZERO) {
  const pts = path.map((p) => Complex.from(p).sub(w));
  let k = 0;
  for (let j = 0; j < pts.length; j++) k += crossing(pts[j], pts[(j + 1) % pts.length]);
  return k;
}

/** +1 when the segment a → b crosses the negative real axis downward (the principal log jumps up by 2πi), −1 upward, else 0. */
function crossing(a, b) {
  if (a.isZero() || b.isZero()) refuse("branch-point", "the path passes through the branch point");
  const up = (c) => c.im.numerator >= 0n; // the axis itself counts as the upper side, as the principal branch has it
  if (up(a) === up(b)) {
    if (a.im.numerator === 0n && b.im.numerator === 0n && a.re.numerator > 0n !== b.re.numerator > 0n) refuse("branch-point", "the path passes through the branch point");
    return 0;
  }
  // where the segment meets the real axis: a + t(b − a) with t = a.im / (a.im − b.im)
  const t = a.im.divide(a.im.subtract(b.im));
  const x = a.re.add(b.re.subtract(a.re).multiply(t));
  if (x.numerator === 0n) refuse("branch-point", "the path passes through the branch point");
  if (x.numerator > 0n) return 0;
  return up(a) ? 1 : -1;
}

/**
 * Carry a branch of log, a root or a power continuously along a polygonal
 * path of exact points, starting from the principal value (or { branch }) at
 * the first point.
 *
 *   continueAlong("log", path)               → { value, turns, values }
 *   continueAlong({ root: 2 }, path)
 *   continueAlong({ pow: "1/3" }, path)
 *
 * `turns` is the number of times the path has gone around 0 (so the value of
 * log at the end is the principal value there plus 2πi·turns), `values` the
 * value at every vertex. `about` moves the branch point from 0 to another exact point.
 */
export function continueAlong(fn, path, opts = {}) {
  const bits = opts.bits ?? 64;
  const about = Complex.from(opts.about ?? 0);
  const pts = path.map((p) => Complex.from(p).sub(about));
  if (!pts.length) throw new RangeError("a path needs at least one point");
  let k = opts.branch ?? 0;
  const at = (c) => {
    const l = addTurns(logRect(Rect.from(c), bits + EXTRA, "above"), k, bits + EXTRA);
    if (fn === "log") return l.round(bits + 4);
    const w = typeof fn === "object" && fn.root ? new Rational(1n, BigInt(fn.root)) : typeof fn === "object" && fn.pow !== undefined ? fn.pow : null;
    if (w === null) throw new TypeError('continueAlong takes "log", { root: n } or { pow: w }');
    return expRect(toRect(w).mul(l), bits + EXTRA).round(bits + 4);
  };
  const values = [at(pts[0])];
  for (let j = 1; j < pts.length; j++) {
    k += crossing(pts[j - 1], pts[j]);
    values.push(at(pts[j]));
  }
  return { value: values[values.length - 1], turns: k - (opts.branch ?? 0), branch: k, values };
}

// ---- as Refinements ---------------------------------------------------------------------------

const FUNCTIONS = { exp, sin, cos, tan, sinh, cosh, tanh, log, sqrt, asin, acos, atan, asinh, acosh, atanh };

/**
 * A function value as a Core Refinement: a float preview at once, then the
 * certified enclosure at `bits`. Works with live() like the rest of RatMath.
 *
 *   refined("exp", "1 + i", { bits: 200 }).quick.value   // [re, im] as doubles
 *   refined("exp", "1 + i", { bits: 200 }).exact()       // a Rect (or a Branched)
 */
export function refined(name, z, opts = {}) {
  const f = FUNCTIONS[name];
  if (!f) throw new RangeError(`no function named ${name}`);
  const bits = opts.bits ?? 64;
  return new Refinement(
    () => {
      const r = f(z, { ...opts, bits: 30 });
      const v = r instanceof Branched ? r.hull : r;
      return { value: v.center.toNumbers(), trust: r instanceof Branched && r.straddles ? "unsure" : "good" };
    },
    function* run() {
      yield { done: 0, of: 1, step: `${name} at ${bits} bits` };
      return { value: f(z, { ...opts, bits }), via: `interval evaluation at ${bits} bits` };
    },
  );
}
