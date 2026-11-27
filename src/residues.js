// Residues, Laurent series, the kinds of isolated singularity, and real
// integrals evaluated by residues with the steps shown.
//
// A formula is first split, by the rules of fractions, into a numerator and a
// denominator that are analytic wherever their own inner parts are. The
// singularities of the formula are then among the zeros of the denominator
// (found and counted by `zeros`), the points where an entire function such as
// exp is fed a pole (essential singularities), and the branch points of
// logarithms and roots. Residues come from proven contour integrals, or for a
// simple pole of P/Q from P/Q′ on the pole's enclosure.

import { Rational } from "@ratmath/core";
import { constants, num, v as variable, op, call, format } from "@ratmath/analytic";
import { Complex, Q, ZERO, ONE, TWO, refuse } from "./complex.js";
import { Rect, Disc, I, iv, sqrtUp, roundUp } from "./enclosure.js";
import { exp as expRect } from "./functions.js";
import { asFunction, parseFormula } from "./taylor.js";
import { Contour, argumentPrinciple, zeros, residueSum, _integrateMany } from "./contour.js";

// ---- formulas as fractions --------------------------------------------------------------------

const ONE_AST = num(1);
const isNum = (n, k) => n.type === "num" && n.value.equals(new Rational(BigInt(k)));
const isOne = (n) => isNum(n, 1);
const same = (a, b) => format(a) === format(b);
const mulA = (a, b) => (isOne(a) ? b : isOne(b) ? a : op("*", a, b));
const powA = (a, k) => (k === 0 || isOne(a) ? ONE_AST : k === 1 ? a : op("^", a, num(k)));
const wholeExponent = (e) => {
  if (e.type === "num" && e.value.denominator === 1n) return Number(e.value.numerator);
  if (e.type === "op" && e.op === "neg" && e.args[0].type === "num" && e.args[0].value.denominator === 1n) return -Number(e.args[0].value.numerator);
  return null;
};
const ENTIRE = new Set(["exp", "sin", "cos", "sinh", "cosh", "rgamma", "zeta1"]);
const BRANCHED = new Set(["log", "ln", "sqrt", "lgamma"]);

/**
 * { n, d } with the formula equal to n/d, and in `info` the inner
 * denominators whose zeros are essential singularities and the inner
 * arguments whose zeros (or poles) are branch points.
 */
function fraction(node, info) {
  switch (node.type) {
    case "num":
    case "const":
    case "var":
      return { n: node, d: ONE_AST };
    case "op": {
      if (node.op === "neg") {
        const a = fraction(node.args[0], info);
        return { n: op("neg", a.n), d: a.d };
      }
      if (node.op === "^") {
        const k = wholeExponent(node.args[1]);
        if (k !== null) {
          const a = fraction(node.args[0], info);
          return k >= 0 ? { n: powA(a.n, k), d: powA(a.d, k) } : { n: powA(a.d, -k), d: powA(a.n, -k) };
        }
        if (node.args[0].type === "const" && node.args[0].name === "e") return fraction(call("exp", node.args[1]), info);
        const [b, e] = [fraction(node.args[0], info), fraction(node.args[1], info)];
        info.branch.push(b.n);
        if (!isOne(b.d)) info.branch.push(b.d);
        if (!isOne(e.d)) info.essential.push(e.d);
        return { n: node, d: ONE_AST };
      }
      const [a, b] = node.args.map((x) => fraction(x, info));
      switch (node.op) {
        case "+":
        case "-":
          if (same(a.d, b.d)) return { n: op(node.op, a.n, b.n), d: a.d };
          return { n: op(node.op, mulA(a.n, b.d), mulA(b.n, a.d)), d: mulA(a.d, b.d) };
        case "*":
          return { n: mulA(a.n, b.n), d: mulA(a.d, b.d) };
        case "/":
          return { n: mulA(a.n, b.d), d: mulA(a.d, b.n) };
      }
      return refuse("unsupported", `operator ${node.op}`);
    }
    case "call": {
      const u = node.args[0];
      if (node.name === "tan") return fraction(op("/", call("sin", u), call("cos", u)), info);
      if (node.name === "tanh") return fraction(op("/", call("sinh", u), call("cosh", u)), info);
      const a = fraction(u, info);
      if (node.name === "gamma" || node.name === "zeta") {
        // Γ = 1/(1/Γ) and ζ(u) = [(u − 1)ζ(u)]/(u − 1), with 1/Γ and (u − 1)ζ(u) entire
        if (!isOne(a.d)) info.essential.push(a.d);
        return node.name === "gamma" ? { n: ONE_AST, d: call("rgamma", u) } : { n: call("zeta1", u), d: op("-", u, num(1)) };
      }
      if (ENTIRE.has(node.name)) {
        if (!isOne(a.d)) info.essential.push(a.d);
        return { n: node, d: ONE_AST };
      }
      if (BRANCHED.has(node.name)) {
        info.branch.push(a.n);
        if (!isOne(a.d)) info.branch.push(a.d);
        return { n: node, d: ONE_AST };
      }
      return refuse("unsupported", `${node.name} is not available here; the functions are exp, log, sqrt, sin, cos, tan, sinh, cosh, tanh, gamma, zeta`);
    }
  }
  return refuse("unsupported", `cannot read a ${node.type} node`);
}
const astOf = (f, what) => {
  if (typeof f === "string") return parseFormula(f);
  if (f && typeof f.type === "string") return f;
  return refuse("bad-argument", `${what} needs the formula as text: a JavaScript function cannot be taken apart into numerator and denominator`);
};
const rename = (node, from, to) => {
  if (node.type === "var") return node.name === from ? variable(to) : node;
  if (node.args) return { ...node, args: node.args.map((a) => rename(a, from, to)) };
  return node;
};
const names = (node, out = new Set()) => {
  if (node.type === "var") out.add(node.name);
  if (node.args) node.args.forEach((a) => names(a, out));
  return out;
};

// ---- polynomials with exact complex coefficients ---------------------------------------------

const C0 = Complex.ZERO;
const ptrim = (a) => {
  while (a.length && a[a.length - 1].isZero()) a.pop();
  return a;
};
const padd = (a, b, sign = 1) => ptrim(Array.from({ length: Math.max(a.length, b.length) }, (_, k) => (a[k] ?? C0)[sign > 0 ? "add" : "sub"](b[k] ?? C0)));
function pmul(a, b) {
  if (!a.length || !b.length) return [];
  const out = new Array(a.length + b.length - 1).fill(C0);
  a.forEach((x, i) => b.forEach((y, j) => (out[i + j] = out[i + j].add(x.mul(y)))));
  return ptrim(out);
}
function pdivmod(a, b) {
  const r = a.slice();
  const q = new Array(Math.max(0, a.length - b.length + 1)).fill(C0);
  const lead = b[b.length - 1];
  for (let k = r.length - b.length; k >= 0; k--) {
    const c = r[k + b.length - 1].div(lead);
    q[k] = c;
    for (let j = 0; j < b.length; j++) r[k + j] = r[k + j].sub(c.mul(b[j]));
  }
  return { q: ptrim(q), r: ptrim(r) };
}
function pgcd(a, b) {
  while (b.length) [a, b] = [b, pdivmod(a, b).r];
  return a.map((c) => c.div(a[a.length - 1]));
}
const pderiv = (a) => ptrim(a.slice(1).map((c, k) => c.scale(new Rational(BigInt(k + 1)))));
/** Horner's rule on anything with mul and add: a Rect, or a Taylor model. */
const pevalOn = (a, z) => {
  let s = null;
  for (let k = a.length - 1; k >= 0; k--) s = s === null ? z.mul(0).add(a[k]) : s.mul(z).add(a[k]);
  return s ?? z.mul(0);
};
const pevalRect = (a, r) => a.reduceRight((s, c) => s.mul(r).add(c), Rect.from(C0));
function ptext(a, x = "z") {
  if (!a.length) return "0";
  const parts = [];
  for (let k = a.length - 1; k >= 0; k--) {
    const c = a[k];
    if (c.isZero()) continue;
    const xs = k === 0 ? "" : k === 1 ? x : `${x}^${k}`;
    let cs = c.toString();
    const compound = !c.isReal() && !c.re.equals(ZERO);
    if (compound) cs = `(${cs})`;
    if (k > 0 && cs === "1") cs = "";
    else if (k > 0 && cs === "-1") cs = "-";
    parts.push(`${cs}${xs}`);
  }
  return parts.join(" + ").replace(/\+ -/g, "- ");
}
/** The polynomial a formula is, in the variable z, or null when it is not one with exact coefficients. */
function polyOf(node) {
  switch (node.type) {
    case "num":
      return ptrim([new Complex(node.value, ZERO)]);
    case "var":
      return node.name === "z" ? [C0, Complex.ONE] : node.name === "i" ? [Complex.I] : null;
    case "op": {
      if (node.op === "^") {
        const k = wholeExponent(node.args[1]);
        const a = polyOf(node.args[0]);
        if (k === null || k < 0 || !a) return null;
        let out = [Complex.ONE];
        for (let j = 0; j < k; j++) out = pmul(out, a);
        return out;
      }
      const args = node.args.map(polyOf);
      if (args.some((a) => a === null)) return null;
      switch (node.op) {
        case "neg": return args[0].map((c) => c.neg());
        case "+": return padd(args[0], args[1]);
        case "-": return padd(args[0], args[1], -1);
        case "*": return pmul(args[0], args[1]);
        case "/": return args[1].length === 1 ? args[0].map((c) => c.div(args[1][0])) : null;
      }
    }
  }
  return null;
}

// ---- 2πi ------------------------------------------------------------------------------------

const twoPi = (bits) => {
  const pi = constants.pi(bits);
  return I.scale(iv(pi.low, pi.high), TWO);
};
/** v / 2πi. */
const over2PiI = (v, bits) => new Rect(I.div(v.im, twoPi(bits)), I.neg(I.div(v.re, twoPi(bits))));
/** 2πi · v. */
const times2PiI = (v, bits) => new Rect(I.neg(I.mul(v.im, twoPi(bits))), I.mul(v.re, twoPi(bits)));
const grow = (r, by) => new Rect(iv(r.re.lo.subtract(by), r.re.hi.add(by)), iv(r.im.lo.subtract(by), r.im.hi.add(by)));

// ---- Laurent series --------------------------------------------------------------------------

/**
 * Laurent coefficients of f about `center`, for the annulus that contains
 * the circle of the given radius: aₙ = (1/2πi) ∮ f(z)/(z − c)^(n+1) dz, each
 * as a proven rectangle. Options: from and to (default −6 … 6), and those of
 * contourIntegral.
 *
 *   laurent("exp(z)/z^2", 0, 1).coefficient(-2)    // 1
 */
export function laurent(f, center, radius, opts = {}) {
  const F = asFunction(f);
  const c = Complex.from(center);
  const r = Q(radius);
  const [from, to] = [opts.from ?? -6, opts.to ?? 6];
  if (!Number.isInteger(from) || !Number.isInteger(to) || from > to || to - from > 200) throw new RangeError("from and to are whole numbers with from ≤ to, at most 200 apart");
  const bits = opts.bits ?? 96;
  const G = (z) => {
    const w = z.sub(c);
    const out = [];
    // f·w^(−n−1) for n = to down to from: start at the lowest power and multiply by w
    let cur = F(z);
    cur = to + 1 >= 0 ? cur.mul(w.inv().pow(to + 1)) : cur.mul(w.pow(-(to + 1)));
    for (let n = to; n >= from; n--) {
      out.push(cur);
      if (n > from) cur = cur.mul(w);
    }
    return out.reverse();
  };
  const values = _integrateMany(G, Contour.circle(c, r), opts).map((v) => over2PiI(v, bits).round(bits));
  return Object.freeze({
    center: c,
    radius: r,
    from,
    to,
    coefficients: Object.freeze(values),
    /** aₙ, for n between from and to. */
    coefficient(n) {
      if (!Number.isInteger(n) || n < from || n > to) throw new RangeError(`the coefficients computed are a_${from} … a_${to}`);
      return values[n - from];
    },
    /** The coefficients of negative index, as [{ n, value }], those certainly not 0 marked. */
    get principalPart() {
      return values.map((value, k) => ({ n: k + from, value, nonzero: !value.containsZero() })).filter((t) => t.n < 0);
    },
    toString(digits = 8) {
      return values.map((v, k) => `a_${k + from} = ${v.toString(digits)}`).join("\n");
    },
  });
}

// ---- residues ---------------------------------------------------------------------------------

/**
 * The residue of f at a point: the integral round the circle of the given
 * radius about it (default 1/4), divided by 2πi. The point is exact, or an
 * enclosure such as `zeros` returns.
 *
 * For a formula given as text the circle is checked: every zero of the
 * formula's denominator inside it must lie within radius/65536 of the point,
 * or the call is refused and a smaller radius is needed. For a JavaScript
 * function nothing can be checked, and the result is the sum of the residues
 * inside the circle.
 */
export function residue(f, point, opts = {}) {
  const at = point instanceof Rect || point instanceof Disc ? Disc.from(point) : new Disc(Complex.from(point));
  const r = Q(opts.radius ?? new Rational(1n, 4n));
  if (r.numerator <= 0n) throw new RangeError("the radius must be positive");
  const big = Contour.circle(at.center, r);
  if (typeof f === "string" || (f && typeof f.type === "string")) {
    const { d } = fraction(astOf(f), { essential: [], branch: [] });
    if (!isOne(d)) {
      let small = r.divide(new Rational(65536n));
      if (small.compareTo(at.radius.multiply(TWO)) < 0) small = at.radius.multiply(TWO);
      if (small.compareTo(r) >= 0) refuse("bad-argument", "the radius is smaller than the enclosure of the point");
      const [inBig, inSmall] = [argumentPrinciple(d, big, opts), argumentPrinciple(d, Contour.circle(at.center, small), opts)];
      if (inBig !== inSmall) refuse("bad-argument", `another singularity lies within ${r} of the point: use a smaller radius`);
    }
  }
  return residueSum(f, big, opts);
}

// ---- singularities ----------------------------------------------------------------------------

/**
 * The singularities of a formula in a rectangle, with their kinds.
 *
 *   singularities("exp(z)/(z^2 + 1)^2", rect([-2, 2], [-2, 2]))
 *   → { singularities: [{ enclosure, kind, order, note }], complete, unresolved }
 *
 * kind is "pole" (order: its order, or the total order of the poles in the
 * enclosure when they could not be separated), "removable" (numerator and
 * denominator vanish together; `order` is then the order of the zero left
 * over, 0 for none), "essential" (an entire function of something with a
 * pole there, as in exp(1/z)) or "branch point" (of a logarithm, root or
 * non-integer power: not isolated, a cut leaves it).
 *
 * What is proved: the numbers of zeros of numerator and denominator inside
 * each enclosure, hence the net order there. That a numerator zero and a
 * denominator zero inside the same tiny enclosure are the same point is not
 * something arithmetic on enclosures can prove; `note` says so when it matters.
 */
export function singularities(f, region, opts = {}) {
  const ast = astOf(f, "singularities");
  const R = Rect.from(region);
  const info = { essential: [], branch: [] };
  const { n, d } = fraction(ast, info);
  const out = [];
  const unresolved = [];
  const margin = new Rational(1n, 1n << 24n);
  const explained = [];
  const scan = (node, kind, note) => {
    const z = zeros(node, R, opts);
    for (const e of z.zeros) {
      out.push({ enclosure: e.enclosure, kind, order: null, note });
      explained.push(grow(e.enclosure, margin));
    }
    return z.unresolved;
  };
  const pending = [];
  for (const e of info.essential) pending.push(...scan(e, "essential", "an entire function of something with a pole here"));
  for (const b of info.branch) pending.push(...scan(b, "branch point", "a logarithm, root or power is taken of something that vanishes or has a pole here; a branch cut leaves this point"));
  if (!isOne(d)) {
    const zd = zeros(d, R, opts);
    pending.push(...zd.unresolved);
    for (const e of zd.zeros) {
      const box = grow(e.enclosure, margin);
      let entry;
      try {
        const frame = Contour.rectangle(box);
        // on so small a frame a function with a zero of order m inside is of size margin^m: the models need bits to match
        const fine = { ...opts, bits: (opts.bits ?? 96) + 26 * (e.multiplicity + 4), maxPieces: 4000 };
        const kD = argumentPrinciple(d, frame, fine);
        let kN = 0;
        let why = "";
        try {
          kN = argumentPrinciple(n, frame, fine);
        } catch (err) {
          if (!err.code) throw err;
          kN = null;
          why = err.message;
        }
        if (kN === null) entry = { enclosure: box, kind: "undetermined", order: null, note: `the denominator vanishes here (order ${kD}), but the numerator could not be bounded: ${why}` };
        else if (kN === 0) entry = { enclosure: e.enclosure, kind: "pole", order: kD, note: e.simple || kD === 1 ? "" : `poles of total order ${kD} in this box; a single pole of that order unless several lie closer together than the tolerance` };
        else if (kN >= kD) entry = { enclosure: box, kind: "removable", order: kN - kD, note: `numerator (${kN}) and denominator (${kD}) both vanish in this box; taken as the same point` };
        else entry = { enclosure: box, kind: "pole", order: kD - kN, note: `the denominator vanishes to order ${kD} and the numerator to order ${kN} in this box; taken as the same point` };
      } catch (err) {
        if (!err.code) throw err;
        entry = { enclosure: box, kind: "undetermined", order: null, note: err.message };
      }
      if (explained.some((x) => x.intersects(entry.enclosure))) continue; // already listed as essential or a branch point
      out.push({ ...entry, boundary: e.boundary });
    }
  }
  for (const u of pending) if (!explained.some((x) => x.intersects(u.region))) unresolved.push(u);
  const key = (s) => s.enclosure.center;
  out.sort((a, b) => key(a).re.compareTo(key(b).re) || key(a).im.compareTo(key(b).im));
  return { singularities: out, complete: unresolved.length === 0, unresolved, numerator: format(n), denominator: format(d) };
}

// ---- real integrals by residues ----------------------------------------------------------------

/** All the roots of a polynomial with their enclosures, proved complete: a box from Cauchy's bound, then `zeros`. */
function allRoots(Qp, opts) {
  const lead = Qp[Qp.length - 1];
  let B = ZERO;
  for (const c of Qp.slice(0, -1)) {
    const m = sqrtUp(c.div(lead).norm());
    if (m.compareTo(B) > 0) B = m;
  }
  B = roundUp(B.add(ONE), 12);
  const found = zeros((z) => pevalOn(Qp, z), new Rect(iv(B.negate(), B), iv(B.negate(), B)), { tolerance: opts.tolerance, bits: opts.bits });
  if (!found.complete || found.count !== Qp.length - 1) refuse("precision", "the poles could not all be isolated");
  return found.zeros;
}

/**
 * Σ of the residues of P(z)/Q(z)·e^(iaz) over the poles that `select` takes.
 * select(enclosure) is true, false, or null when the enclosure does not decide.
 */
function sumResidues(P, Qp, a, select, onEdge, opts) {
  const bits = opts.bits ?? 64;
  const roots = allRoots(Qp, opts);
  const dQ = pderiv(Qp);
  const withExp = (v, z) => (a ? v.mul(expRect(z.timesI().scale(a), bits + 16)) : v);
  const poles = [];
  let sum = Rect.from(C0);
  for (const root of roots) {
    const take = select(root.enclosure);
    if (take === null) refuse("pole", onEdge(root.enclosure));
    let res = null;
    if (take && root.simple) {
      const E = root.enclosure;
      res = withExp(pevalRect(P, E), E).div(pevalRect(dQ, E));
    } else if (take) {
      // a multiple pole: a small circle round the box, clear of the other poles
      const c = root.enclosure.center;
      const [x, y] = c.toNumbers();
      let rho = 1;
      for (const o of roots) {
        if (o === root) continue;
        const [u, w] = o.enclosure.center.toNumbers();
        rho = Math.min(rho, Math.hypot(u - x, w - y) / 2);
      }
      const radius = new Rational(BigInt(Math.floor(rho * 2 ** 30)), 1n << 30n);
      const ia = a ? new Complex(ZERO, a) : null;
      const G = (z) => {
        const q = pevalOn(P, z).div(pevalOn(Qp, z));
        return ia ? q.mul(expRect(z.mul(ia))) : q;
      };
      res = residueSum(G, Contour.circle(new Complex(roundUp(c.re, 40), roundUp(c.im, 40)), radius), { order: 16, bits: 128, tolerance: new Rational(1n, 1n << BigInt(bits)) });
    }
    if (res) sum = sum.add(res);
    poles.push({ enclosure: root.enclosure, order: root.multiplicity, taken: Boolean(take), residue: res && res.round(bits + 8) });
  }
  return { poles, sum };
}
const realCoefficients = (p) => p.every((c) => c.isReal());
const rationalFunction = (ast, what) => {
  const { n, d } = fraction(ast, { essential: [], branch: [] });
  let [P, Qp] = [polyOf(n), polyOf(d)];
  if (!P || !Qp || !Qp.length) refuse("unsupported", `${what} needs a quotient of polynomials with exact coefficients`);
  const g = pgcd(P.slice(), Qp.slice());
  if (g.length > 1) [P, Qp] = [pdivmod(P, g).q, pdivmod(Qp, g).q];
  return { P, Qp };
};
const show = (r, digits = 12) => r.toString(digits);
const poleLines = (poles) => poles.filter((p) => p.taken).map((p) => `at ${show(p.enclosure, 10)}${p.order > 1 ? ` (order ${p.order})` : ""}: residue ${show(p.residue, 10)}`);

/**
 * ∫ f(x) dx over the whole real line by residues, for f = P/Q a quotient of
 * polynomials with deg Q ≥ deg P + 2 and no pole on the real axis: 2πi times
 * the sum of the residues in the upper half-plane.
 *
 * With { fourier: a } the integrand is f(x)·e^(iax), and deg Q ≥ deg P + 1 is
 * enough (Jordan's lemma); for real f the result's `cos` and `sin` are then
 * ∫ f(x) cos(ax) dx and ∫ f(x) sin(ax) dx.
 *
 * Returns { value (a Rect), real (the real part's interval, when f is real),
 * cos, sin, poles, steps }. `steps` is the argument, a list of { title, detail }.
 */
export function realLineIntegral(f, opts = {}) {
  let ast = astOf(f, "realLineIntegral");
  const used = names(ast);
  const x = opts.variable ?? (used.has("x") ? "x" : "z");
  if (x !== "z") {
    if (used.has("z")) refuse("bad-argument", `the formula uses both ${x} and z`);
    ast = rename(ast, x, "z");
  }
  const { P, Qp } = rationalFunction(ast, "realLineIntegral");
  const a = opts.fourier === undefined ? null : Q(opts.fourier);
  const fourier = a !== null && a.numerator !== 0n;
  const ia = fourier ? (a.equals(ONE) ? "i" : a.equals(ONE.negate()) ? "-i" : `${a}i`) : "";
  const bits = opts.bits ?? 64;
  const steps = [];
  const [dp, dq] = [P.length - 1, Qp.length - 1];
  steps.push({ title: "The integrand", detail: `f(${x}) = (${ptext(P, x)}) / (${ptext(Qp, x)})${fourier ? `, times e^(${ia}${x})` : ""}, in lowest terms` });
  if (!P.length) return Object.freeze({ value: Rect.from(C0), real: iv(ZERO, ZERO), poles: [], steps });
  const need = fourier ? 1 : 2;
  if (dq - dp < need) refuse("bad-argument", `the integral does not converge: the denominator's degree (${dq}) must exceed the numerator's (${dp}) by at least ${need}`);
  steps.push({
    title: "Why the arc can be dropped",
    detail: fourier
      ? `deg Q − deg P = ${dq - dp} ≥ 1, so f → 0 at infinity and by Jordan's lemma the integral of f(z)e^(${ia}z) over the large semicircle in the ${a.numerator > 0n ? "upper" : "lower"} half-plane tends to 0${dq - dp === 1 ? " (the real integral converges only as a limit of symmetric ranges)" : ""}`
      : `deg Q − deg P = ${dq - dp} ≥ 2, so |f(z)| ≤ C/|z|² and the integral over the semicircle of radius R is at most πC/R, which tends to 0`,
  });
  const upper = !fourier || a.numerator > 0n;
  const select = (E) => (E.im.lo.numerator > 0n ? upper : E.im.hi.numerator < 0n ? !upper : null);
  const { poles, sum } = sumResidues(P, Qp, fourier ? a : null, select, (E) => `a pole lies on the real axis, or too close to it to tell, near ${show(E, 8)}: the integral does not exist as it stands`, opts);
  steps.push({ title: `The poles in the ${upper ? "upper" : "lower"} half-plane`, detail: poleLines(poles).join("; ") || "none" });
  let value = times2PiI(sum, bits + 16);
  if (!upper) value = value.neg();
  value = value.round(bits);
  steps.push({ title: "The residue theorem", detail: `the integral is ${upper ? "" : "−"}2πi × (the sum of those residues) = ${show(value)}` });
  const out = { value, poles, steps };
  if (realCoefficients(P) && realCoefficients(Qp)) {
    if (fourier) Object.assign(out, { cos: value.re, sin: value.im });
    else out.real = value.re;
    steps.push({ title: "Real and imaginary parts", detail: fourier ? `f is real on the real axis, so ∫ f(${x}) cos(${a.equals(ONE) ? "" : a}${x}) d${x} is the real part and ∫ f(${x}) sin(${a.equals(ONE) ? "" : a}${x}) d${x} the imaginary part` : "f is real on the real axis, so the integral is the real part; the imaginary part encloses 0" });
  }
  return Object.freeze(out);
}

/**
 * ∫₀^2π R(cos t, sin t) dt by residues, for R a quotient of polynomials in
 * cos t and sin t (cos 2t and the like are understood): with z = e^(it) the
 * integral becomes one round the unit circle, 2πi times the sum of the
 * residues inside it.
 *
 *   trigIntegral("1/(2 + cos(t))").real      // 2π/√3
 */
export function trigIntegral(f, opts = {}) {
  const ast = astOf(f, "trigIntegral");
  const used = names(ast);
  const t = opts.variable ?? ["t", "theta", "x"].find((name) => used.has(name)) ?? "t";
  const bits = opts.bits ?? 64;
  const Z = variable("z");
  const multiple = (u) => {
    if (u.type === "var" && u.name === t) return 1;
    if (u.type === "op" && u.op === "neg") {
      const k = multiple(u.args[0]);
      return k === null ? null : -k;
    }
    if (u.type === "op" && u.op === "*") {
      const [p, q] = u.args;
      const k = wholeExponent(p) !== null && q.type === "var" && q.name === t ? wholeExponent(p) : wholeExponent(q) !== null && p.type === "var" && p.name === t ? wholeExponent(q) : null;
      return k;
    }
    return null;
  };
  const zk = (k) => op("^", Z, num(k));
  const swap = (node) => {
    if (node.type === "var" && node.name !== "i") return refuse("bad-argument", `${node.name} appears outside cos and sin: the integrand must be built from cos(${t}) and sin(${t})`);
    if (node.type === "call" && ["cos", "sin", "tan"].includes(node.name)) {
      const k = multiple(node.args[0]);
      if (k === null) return refuse("bad-argument", `${node.name} must be of ${t} or a whole multiple of it`);
      const c = op("/", op("+", zk(k), zk(-k)), num(2));
      const s = op("/", op("-", zk(k), zk(-k)), op("*", num(2), variable("i")));
      return node.name === "cos" ? c : node.name === "sin" ? s : op("/", s, c);
    }
    if (node.type === "call") return refuse("unsupported", `${node.name} is not available here: only cos, sin and tan`);
    if (node.args) return { ...node, args: node.args.map(swap) };
    return node;
  };
  const g = op("/", swap(ast), op("*", variable("i"), Z));
  const { P, Qp } = rationalFunction(g, "trigIntegral");
  const steps = [
    { title: "Substitute z = e^(it)", detail: `cos t = (z + 1/z)/2, sin t = (z − 1/z)/(2i), dt = dz/(iz): the integral is ∮ g(z) dz round the unit circle with g(z) = (${ptext(P)}) / (${ptext(Qp)})` },
  ];
  const select = (E) => {
    const n = E.normRange();
    return n.hi.compareTo(ONE) < 0 ? true : n.lo.compareTo(ONE) > 0 ? false : null;
  };
  const { poles, sum } = sumResidues(P, Qp, null, select, (E) => `a pole lies on the unit circle, or too close to it to tell, near ${show(E, 8)}: the integrand is unbounded on [0, 2π]`, opts);
  steps.push({ title: "The poles inside the unit circle", detail: poleLines(poles).join("; ") || "none" });
  const value = times2PiI(sum, bits + 16).round(bits);
  steps.push({ title: "The residue theorem", detail: `the integral is 2πi × (the sum of those residues) = ${show(value)}` });
  return Object.freeze({ value, real: value.re, poles, steps });
}
