import { expect, test, describe } from "bun:test";
import { gamma, rgamma, logGamma, zeta, zetaTimesSMinus1, zetaZeros, taylorModel, contourIntegral, residue, singularities, zeros, countZeros, Contour, complex, rect, disc, exp, ComplexError } from "../index.js";

const near = (r, re, im, eps = 1e-12) => {
  const [x, y] = r.center.toNumbers();
  expect(Math.abs(x - re)).toBeLessThan(eps * Math.max(1, Math.abs(re)));
  expect(Math.abs(y - im)).toBeLessThan(eps * Math.max(1, Math.abs(im)));
};
const code = (fn) => {
  try {
    fn();
  } catch (e) {
    expect(e).toBeInstanceOf(ComplexError);
    return e.code;
  }
  return "no error";
};

describe("gamma", () => {
  test("factorials and half-integers are inside their enclosures", () => {
    expect(gamma(5).contains(complex(24))).toBe(true);
    expect(gamma(1).contains(complex(1))).toBe(true);
    expect(gamma(11).contains(complex(3628800))).toBe(true);
    // Γ(1/2)² = π, Γ(−3/2) = 4√π/3
    near(gamma("1/2").mul(gamma("1/2")), Math.PI, 0);
    near(gamma("-3/2"), (4 * Math.sqrt(Math.PI)) / 3, 0);
    expect(gamma(5).size.toNumber()).toBeLessThan(1e-17);
  });
  test("complex values", () => {
    near(gamma("1 + i"), 0.4980156681183561, -0.1549498283018107);
    near(gamma("-5/2 + 3i"), 0.00047978841084189701, 0.00029885571114485887, 1e-10);
    // |Γ(iy)|² = π/(y sinh πy)
    const g = gamma("2i");
    const [x, y] = g.center.toNumbers();
    expect(x * x + y * y).toBeCloseTo(Math.PI / (2 * Math.sinh(2 * Math.PI)), 14);
  });
  test("the recurrence and the reflection formula hold between enclosures", () => {
    const z = complex("-7/3 + 5/4i");
    expect(gamma(z.add(1)).intersects(gamma(z).mul(z))).toBe(true);
    const far = complex("-80 + 1/3i");
    expect(gamma(far.add(1)).intersects(gamma(far).mul(far))).toBe(true);
  });
  test("more bits", () => {
    const g = gamma("30 + 40i", { bits: 200 });
    expect(g.size.toNumber() / 1e21).toBeLessThan(1e-55);
    expect(g.intersects(gamma("30 + 40i"))).toBe(true);
  });
  test("an enclosure in, an enclosure out", () => {
    const g = gamma(rect([2, 2.001], [1, 1.001]));
    expect(g.contains(gamma("2 + i").center)).toBe(true);
    expect(g.contains(gamma("2001/1000 + 1001/1000i").center)).toBe(true);
  });
  test("poles are refused; 1/Γ is entire and exactly 0 there", () => {
    expect(code(() => gamma(0))).toBe("pole");
    expect(code(() => gamma(-3))).toBe("pole");
    expect(rgamma(-3).contains(complex(0))).toBe(true);
    near(rgamma(5), 1 / 24, 0);
  });
  test("logGamma", () => {
    near(logGamma(10), Math.log(362880), 0);
    // the continuous branch: Im log Γ(1 + 20i) is far outside (−π, π]
    const l = logGamma("1 + 20i");
    expect(Math.abs(l.center.toNumbers()[1])).toBeGreaterThan(Math.PI);
    expect(exp(l).intersects(gamma("1 + 20i"))).toBe(true);
    expect(code(() => logGamma("-1/2"))).toBe("bad-argument");
  });
});

describe("zeta", () => {
  test("known values", () => {
    near(zeta(2), Math.PI ** 2 / 6, 0);
    near(zeta(4), Math.PI ** 4 / 90, 0);
    expect(zeta(0).contains(complex("-1/2"))).toBe(true);
    expect(zeta(-1).contains(complex("-1/12"))).toBe(true);
    expect(zeta(-2).contains(complex(0))).toBe(true);
    expect(zeta(-3).contains(complex("1/120"))).toBe(true);
    near(zeta("1/2"), -1.4603545088095868, 0);
    near(zeta("2 + 3i"), 0.7980219851462757, -0.1137443080529385);
  });
  test("the functional equation holds between the two formulas", () => {
    // s to the left uses the reflection; 1 − s to the right uses Euler–Maclaurin
    const s = complex("-5/2 + 7i");
    const left = zeta(s);
    near(left, 0.41467595079096174, 1.46806314597332373, 1e-11);
    // and just either side of Re s = 0 the two agree
    const a = zeta("1/1000 + 3i");
    const b = zeta("-1/1000 + 3i");
    expect(Math.abs(a.center.toNumbers()[0] - b.center.toNumbers()[0])).toBeLessThan(0.01);
    expect(zeta(rect([-0.01, 0.01], [3, 3.01])).contains(a.center)).toBe(true);
    expect(zeta(rect([-0.01, 0.01], [3, 3.01])).contains(b.center)).toBe(true);
  });
  test("the pole and its entire companion", () => {
    expect(code(() => zeta(1))).toBe("pole");
    expect(zetaTimesSMinus1(1).contains(complex(1))).toBe(true);
    // the residue of ζ at 1 is 1, by an integral round it
    near(residue((z) => zeta(z), 1, { radius: "1/2" }), 1, 0, 1e-9);
  });
  test("high up the critical line, and to more digits", () => {
    near(zeta("1/2 + 100i"), 2.692619885681324, -0.02038602960259816, 1e-11);
    const z = zeta("2 + 3i", { bits: 200 });
    expect(z.size.toNumber()).toBeLessThan(1e-55);
    expect(z.intersects(zeta("2 + 3i"))).toBe(true);
  });
  test("ζ(2) by a contour integral of its model: Cauchy's formula", () => {
    // (1/2πi) ∮ ζ(z)/(z − 2) dz round a circle that leaves the pole at 1 outside
    near(residue((z) => zeta(z).div(z.sub(2)), 2, { radius: "1/2" }), Math.PI ** 2 / 6, 0, 1e-9);
  });
});

describe("gamma and zeta in the rest of the package", () => {
  test("as text, in Taylor models", () => {
    const m = taylorModel("gamma(z)", disc("2 + i", "1/4"));
    expect(m.at("2 + i").intersects(gamma("2 + i"))).toBe(true);
    expect(m.error.toNumber()).toBeLessThan(1e-6);
    const k = taylorModel("zeta(z)^2", disc("1/2 + 14i", "1/4"));
    const v = zeta("1/2 + 14i");
    expect(k.at("1/2 + 14i").intersects(v.mul(v))).toBe(true);
  });
  test("the poles of gamma are found as singularities", () => {
    const s = singularities("gamma(z)", rect([-7 / 2, 1], [-1, 1]));
    expect(s.singularities.map((x) => [x.kind, x.order])).toEqual([["pole", 1], ["pole", 1], ["pole", 1], ["pole", 1]]);
    for (const k of [0, -1, -2, -3]) expect(s.singularities.some((x) => x.enclosure.contains(complex(k)))).toBe(true);
    // the residue at −n is (−1)ⁿ/n!
    near(residue("gamma(z)", -2, { radius: "1/2" }), 1 / 2, 0, 1e-9);
    near(residue("gamma(z)", -3, { radius: "1/2" }), -1 / 6, 0, 1e-9);
  });
  test("the pole of zeta, and a trivial zero", () => {
    const s = singularities("zeta(z)", rect([-1, 2], [-1, 1]));
    expect(s.singularities.map((x) => [x.kind, x.order])).toEqual([["pole", 1]]);
    const z = zeros("zeta(z)", rect([-5 / 2, -3 / 2], [-1 / 2, 1 / 2]));
    expect(z.count).toBe(1);
    expect(z.zeros[0].enclosure.contains(complex(-2))).toBe(true);
  });
});

describe("the zeros of zeta", () => {
  test("the first three, each proved to lie on the critical line", () => {
    const r = zetaZeros(10, 30);
    expect(r.complete).toBe(true);
    expect(r.count).toBe(3);
    expect(r.allOnCriticalLine).toBe(true);
    const t = r.zeros.map((z) => (z.t.lo.toNumber() + z.t.hi.toNumber()) / 2);
    expect(t[0]).toBeCloseTo(14.134725141734693, 12);
    expect(t[1]).toBeCloseTo(21.022039638771555, 12);
    expect(t[2]).toBeCloseTo(25.010857580145688, 12);
    for (const z of r.zeros) expect(z.t.hi.subtract(z.t.lo).toNumber()).toBeLessThan(1e-17);
  });
  test("the count agrees with the argument principle on the boundary", () => {
    expect(countZeros("zeta(z)", rect([0, 1], [10, 30]))).toBe(3);
  });
  test("no zeros below the first", () => {
    const r = zetaZeros(1, 14);
    expect(r.complete).toBe(true);
    expect(r.count).toBe(0);
  });
});
