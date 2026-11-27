// The Cauchy–Riemann equations for f = u + iv given by its real and imaginary
// parts as formulas in x and y. The partial derivatives are symbolic (analytic's
// differentiate, simplify and tidy); a residual that reduces to 0 is zero
// wherever both sides are defined.
// One that does not is evaluated at sample points with certified enclosures: an
// enclosure that excludes 0 is a witness that the equations fail there.

import { Rational } from "@ratmath/core";
import { parse, format, differentiate, simplify, tidy, field } from "@ratmath/analytic";
import { refuse } from "./complex.js";

const SAMPLES = ["1/3", "1/2", "2/3", "1", "3/2", "2", "-1/2", "-1", "-3/2", "5/7"];

const isZero = (ast) => ast.type === "num" && ast.value.numerator === 0n;

function residual(text, [x, y]) {
  // a nonzero value at a sample point, proven by an enclosure that excludes 0
  if (text === "0") return null;
  const F = field(text, [x, y]);
  for (const a of SAMPLES) {
    for (const b of SAMPLES) {
      try {
        const v = F.value([a, b], { digits: 10 }).exact();
        if (v.low.greaterThan(new Rational(0)) || v.high.lessThan(new Rational(0))) return { point: [a, b], value: v };
      } catch {
        // outside the domain at this point: try the next one
      }
    }
  }
  return undefined;
}

/**
 * Check the Cauchy–Riemann equations u_x = v_y and u_y = −v_x for f = u + iv.
 * @param {string} u  real part, a formula in x and y
 * @param {string} v  imaginary part
 * @param {{ variables?: [string, string] }} [options]
 * @returns {{ holds: boolean|null, reason: string, ux: string, uy: string, vx: string, vy: string,
 *   residuals: [string, string], derivative: string|null, witness?: object }}
 *   `holds` is true when both residuals simplify to 0, false when one is proven
 *   nonzero at a point (the witness), null when neither is shown. When they hold,
 *   `derivative` is f′ = u_x + i·v_x, and `derivativeParts` its real and
 *   imaginary parts as formulas.
 */
export function cauchyRiemann(u, v, { variables = ["x", "y"] } = {}) {
  if (!Array.isArray(variables) || variables.length !== 2) refuse("bad-argument", "give two variable names, such as [\"x\", \"y\"]");
  const [x, y] = variables;
  const U = parse(String(u));
  const V = parse(String(v));
  const d = (ast, name) => simplify(differentiate(ast, name));
  const ux = d(U, x), uy = d(U, y), vx = d(V, x), vy = d(V, y);
  // tidy: a rational function of x, y and the calls in lowest terms, so log/atan pairs cancel
  const normal = (text) => {
    const s = simplify(parse(text));
    try {
      return tidy(s);
    } catch {
      return s;
    }
  };
  const r1 = normal(`(${format(ux)}) - (${format(vy)})`);
  const r2 = normal(`(${format(uy)}) + (${format(vx)})`);
  const out = {
    ux: format(ux), uy: format(uy), vx: format(vx), vy: format(vy),
    residuals: [format(r1), format(r2)],
  };
  if (isZero(r1) && isZero(r2)) {
    return { holds: true, reason: `u_${x} = v_${y} = ${out.ux} and u_${y} = −v_${x} = ${out.uy} identically`, ...out, derivative: `${format(normal(out.ux))} + i·(${format(normal(out.vx))})`, derivativeParts: [format(normal(out.ux)), format(normal(out.vx))] };
  }
  for (const [i, r] of [[0, r1], [1, r2]]) {
    if (isZero(r)) continue;
    const w = residual(format(r), variables);
    if (w) {
      const eq = i === 0 ? `u_${x} − v_${y}` : `u_${y} + v_${x}`;
      return { holds: false, reason: `${eq} = ${format(r)}, which is not 0 at (${x}, ${y}) = (${w.point.join(", ")})`, ...out, derivative: null, witness: { equation: eq, point: w.point, value: w.value } };
    }
  }
  return { holds: null, reason: `the residuals ${out.residuals.join(" and ")} do not simplify to 0, and no sample point showed them nonzero`, ...out, derivative: null };
}
