import { describe, expect, test } from "bun:test";
import { field } from "@ratmath/analytic";
import { cauchyRiemann, complex } from "../index.js";

describe("cauchyRiemann", () => {
  test("z² = (x² − y²) + 2xy i", () => {
    const r = cauchyRiemann("x^2 - y^2", "2*x*y");
    expect(r.holds).toBe(true);
    expect([r.ux, r.vy, r.uy, r.vx]).toEqual(["2*x", "2*x", "-2*y", "2*y"]);
    expect(r.derivative).toBe("2*x + i·(2*y)");
  });

  test("entire functions hold; transcendental and rational", () => {
    expect(cauchyRiemann("exp(x)*cos(y)", "exp(x)*sin(y)").holds).toBe(true);
    expect(cauchyRiemann("sin(x)*cosh(y)", "cos(x)*sinh(y)").holds).toBe(true);
    expect(cauchyRiemann("x^3 - 3*x*y^2", "3*x^2*y - y^3").holds).toBe(true);
    // 1/z and the principal log away from the cut
    expect(cauchyRiemann("x/(x^2+y^2)", "-y/(x^2+y^2)").holds).toBe(true);
    expect(cauchyRiemann("log(x^2+y^2)/2", "atan(y/x)").holds).toBe(true);
  });

  test("conj(z), |z|² and 1/conj(z) fail, with a witness point", () => {
    for (const [u, v] of [["x", "-y"], ["x^2+y^2", "0"], ["x/(x^2+y^2)", "y/(x^2+y^2)"]]) {
      const r = cauchyRiemann(u, v);
      expect(r.holds).toBe(false);
      expect(r.witness.point.length).toBe(2);
      // the witness enclosure excludes 0
      const { low, high } = r.witness.value;
      expect(low.numerator > 0n || high.numerator < 0n).toBe(true);
    }
  });

  test("other variable names", () => {
    expect(cauchyRiemann("a^2 - b^2", "2*a*b", { variables: ["a", "b"] }).holds).toBe(true);
    expect(() => cauchyRiemann("x", "y", { variables: ["x"] })).toThrow();
  });

  test("the derivative agrees with the complex function's own derivative (independent check)", () => {
    // f = z³: f′ = 3z²; at z = 1/2 + i/3 that is 3(1/4 − 1/9) + 3·(1/3) i = 5/12 + i
    const r = cauchyRiemann("x^3 - 3*x*y^2", "3*x^2*y - y^3");
    const [re, im] = r.derivativeParts.map((t) => field(t, ["x", "y"]).value(["1/2", "1/3"]).exact());
    const z = complex("1/2", "1/3");
    const d = z.mul(z).mul(complex(3, 0));
    expect(re.low.equals(d.re) && re.high.equals(d.re)).toBe(true);
    expect(im.low.equals(d.im) && im.high.equals(d.im)).toBe(true);
  });
});
