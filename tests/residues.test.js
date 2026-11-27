import { expect, test, describe } from "bun:test";
import { laurent, residue, singularities, realLineIntegral, trigIntegral, rect, complex, zeros, ComplexError } from "../index.js";

const near = (r, re, im, eps = 1e-9) => {
  const [x, y] = r.center.toNumbers();
  expect(Math.abs(x - re)).toBeLessThan(eps);
  expect(Math.abs(y - im)).toBeLessThan(eps);
};
const mid = (i) => (i.lo.toNumber() + i.hi.toNumber()) / 2;
const code = (fn) => {
  try {
    fn();
  } catch (e) {
    expect(e).toBeInstanceOf(ComplexError);
    return e.code;
  }
  return "no error";
};

describe("Laurent series", () => {
  test("exp(z)/z² about 0", () => {
    const L = laurent("exp(z)/z^2", 0, 1, { from: -3, to: 3 });
    const want = { "-3": 0, "-2": 1, "-1": 1, 0: 1 / 2, 1: 1 / 6, 2: 1 / 24, 3: 1 / 120 };
    for (let n = -3; n <= 3; n++) near(L.coefficient(n), want[n], 0);
    expect(L.principalPart.filter((t) => t.nonzero).map((t) => t.n)).toEqual([-2, -1]);
    expect(L.coefficient(-3).containsZero()).toBe(true);
    expect(() => L.coefficient(9)).toThrow(RangeError);
  });
  test("an essential singularity has an endless principal part", () => {
    const L = laurent("exp(1/z)", 0, 1, { from: -5, to: 0 });
    let f = 1;
    for (let n = 0; n <= 5; n++) {
      near(L.coefficient(-n), 1 / f, 0);
      f *= n + 1;
    }
    expect(L.principalPart.every((t) => t.nonzero)).toBe(true);
  });
  test("the series depends on the annulus", () => {
    const f = "1/((z - 1)*(z - 2))";
    // inside |z| < 1: Σ (1 − 2^−(n+1)) zⁿ; between 1 and 2: −Σ_{n<0} zⁿ − Σ_{n≥0} zⁿ/2^(n+1)
    const inner = laurent(f, 0, "1/2", { from: -2, to: 2 });
    near(inner.coefficient(0), 1 / 2, 0);
    near(inner.coefficient(1), 3 / 4, 0);
    expect(inner.coefficient(-1).containsZero()).toBe(true);
    const middle = laurent(f, 0, "3/2", { from: -2, to: 2 });
    near(middle.coefficient(-1), -1, 0);
    near(middle.coefficient(-2), -1, 0);
    near(middle.coefficient(0), -1 / 2, 0);
    near(middle.coefficient(1), -1 / 4, 0);
  });
});

describe("residues", () => {
  test("at poles of several orders", () => {
    near(residue("exp(z)/z^3", 0), 1 / 2, 0);
    near(residue("1/(z^2 + 1)", "i"), 0, -1 / 2);
    near(residue("cos(z)/sin(z)", 0, { radius: 1 }), 1, 0);
    near(residue("exp(1/z)", 0), 1, 0);
  });
  test("at an enclosure from zeros", () => {
    const pole = zeros("z^2 - 2", rect([1, 2], [-1, 1])).zeros[0];
    // 1/(z² − 2) at √2: 1/(2√2)
    near(residue("1/(z^2 - 2)", pole.enclosure, { radius: "1/2" }), 1 / (2 * Math.SQRT2), 0);
  });
  test("a circle that takes in a second singularity is refused", () => {
    expect(code(() => residue("1/(z*(z - 1))", 0, { radius: 2 }))).toBe("bad-argument");
    near(residue("1/(z*(z - 1))", 0, { radius: "1/2" }), -1, 0);
  });
  test("a JavaScript function is integrated as it is", () => {
    near(residue((z) => z.inv(), 0), 1, 0);
  });
});

describe("singularities", () => {
  const kinds = (f, R) => singularities(f, R).singularities.map((s) => [s.kind, s.order]);
  test("poles with their orders", () => {
    const s = singularities("exp(z)/(z^2 + 1)^2", rect([-2, 2], [-2, 2]));
    expect(s.complete).toBe(true);
    expect(s.singularities.map((x) => [x.kind, x.order])).toEqual([["pole", 2], ["pole", 2]]);
    expect(s.singularities.some((x) => x.enclosure.contains(complex("i")))).toBe(true);
    expect(kinds("tan(z)", rect([-4, 4], [-1, 1]))).toEqual([["pole", 1], ["pole", 1]]);
    expect(kinds("1/(exp(z) - 1)", rect([-1, 1], [-7, 7]))).toEqual([["pole", 1], ["pole", 1], ["pole", 1]]);
  });
  test("removable singularities and cancelled orders", () => {
    expect(kinds("sin(z)/z", rect([-1, 1], [-1, 1]))).toEqual([["removable", 0]]);
    expect(kinds("sin(z)^3/z", rect([-1, 1], [-1, 1]))).toEqual([["removable", 2]]);
    expect(kinds("(1 - cos(z))/z^5", rect([-1, 1], [-1, 1]))).toEqual([["pole", 3]]);
  });
  test("essential singularities and branch points", () => {
    expect(kinds("exp(1/z)", rect([-1, 1], [-1, 1]))).toEqual([["essential", null]]);
    expect(kinds("sin(1/(z - 1/2)) + 1/z", rect([-1, 1], [-1, 1]))).toEqual([["pole", 1], ["essential", null]]);
    const s = singularities("sqrt(z - 1/2)", rect([-1, 1], [-1, 1]));
    expect(s.singularities.map((x) => x.kind)).toEqual(["branch point"]);
  });
  test("an entire function has none", () => {
    const s = singularities("exp(z) + z^2", rect([-1, 1], [-1, 1]));
    expect(s.singularities).toEqual([]);
    expect(s.complete).toBe(true);
  });
  test("a JavaScript function cannot be taken apart", () => {
    expect(code(() => singularities((z) => z.inv(), rect([-1, 1], [-1, 1])))).toBe("bad-argument");
  });
});

describe("real integrals by residues", () => {
  test("∫ dx/(x² + 1) = π, with the steps", () => {
    const r = realLineIntegral("1/(x^2 + 1)");
    expect(mid(r.real)).toBeCloseTo(Math.PI, 14);
    expect(r.real.hi.subtract(r.real.lo).toNumber()).toBeLessThan(1e-17);
    expect(r.value.im.lo.numerator <= 0n && r.value.im.hi.numerator >= 0n).toBe(true);
    expect(r.steps.map((s) => s.title)).toEqual(["The integrand", "Why the arc can be dropped", "The poles in the upper half-plane", "The residue theorem", "Real and imaginary parts"]);
    expect(r.poles.filter((p) => p.taken).length).toBe(1);
  });
  test("textbook values", () => {
    expect(mid(realLineIntegral("1/(x^2 + 1)^2").real)).toBeCloseTo(Math.PI / 2, 13);
    expect(mid(realLineIntegral("x^2/(x^4 + 1)").real)).toBeCloseTo(Math.PI / Math.SQRT2, 13);
    expect(mid(realLineIntegral("1/(x^4 + 5x^2 + 4)").real)).toBeCloseTo(Math.PI / 6, 13);
    expect(mid(realLineIntegral("1/(z^2 + 2z + 2)").real)).toBeCloseTo(Math.PI, 13);
  });
  test("common factors are cancelled first", () => {
    expect(mid(realLineIntegral("(x^2 - 1)/((x^2 - 1)*(x^2 + 4))").real)).toBeCloseTo(Math.PI / 2, 13);
  });
  test("Fourier integrals by Jordan's lemma", () => {
    // ∫ cos x/(x² + 1) = π/e;  ∫ x sin 2x/(x² + 1) = π e^−2
    const a = realLineIntegral("1/(x^2 + 1)", { fourier: 1 });
    expect(mid(a.cos)).toBeCloseTo(Math.PI / Math.E, 13);
    expect(mid(a.sin)).toBeCloseTo(0, 13);
    const b = realLineIntegral("x/(x^2 + 1)", { fourier: 2 });
    expect(mid(b.sin)).toBeCloseTo(Math.PI * Math.exp(-2), 13);
    // a negative frequency closes in the lower half-plane
    const c = realLineIntegral("1/(x^2 + 1)", { fourier: -1 });
    expect(mid(c.cos)).toBeCloseTo(Math.PI / Math.E, 13);
  });
  test("what is refused", () => {
    expect(code(() => realLineIntegral("1/(x^2 - 1)"))).toBe("pole");
    expect(code(() => realLineIntegral("x/(x^2 + 1)"))).toBe("bad-argument");
    expect(code(() => realLineIntegral("exp(x)/(x^2 + 1)"))).toBe("unsupported");
  });
  test("trigonometric integrals over a period", () => {
    expect(mid(trigIntegral("1/(2 + cos(t))").real)).toBeCloseTo((2 * Math.PI) / Math.sqrt(3), 13);
    expect(mid(trigIntegral("cos(t)^4").real)).toBeCloseTo((3 * Math.PI) / 4, 13);
    expect(mid(trigIntegral("cos(2t)/(5 - 4cos(t))").real)).toBeCloseTo(Math.PI / 6, 13);
    expect(mid(trigIntegral("1/(1 + sin(theta)^2)").real)).toBeCloseTo(Math.PI * Math.SQRT2, 13);
    expect(mid(trigIntegral("sin(t)").real)).toBeCloseTo(0, 13);
    expect(code(() => trigIntegral("1/cos(t)"))).toBe("pole");
    expect(code(() => trigIntegral("t*cos(t)"))).toBe("bad-argument");
  });
});
