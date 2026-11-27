import { expect, test, describe } from "bun:test";
import { Rational } from "@ratmath/core";
import { taylorModel, CTM, complex, disc, rect, exp, sin, cos, tan, sinh, cosh, tanh, log, sqrt, root, pow, Rect, ComplexError } from "../index.js";

// a few exact points of a disc, the centre and the rim included
const pointsOf = (D) => {
  const out = [D.center, D.center.add(complex(D.radius)), D.center.sub(complex(0, D.radius))];
  for (const [a, b] of [[3, 5], [-4, 7], [4, -9], [-2, -3], [1, 8]]) {
    out.push(D.center.add(complex(new Rational(BigInt(a), 10n), new Rational(BigInt(b), 10n)).scale(D.radius)));
  }
  return out;
};
const meets = (a, b) => Rect.from(a).intersects(Rect.from(b));

describe("Taylor models enclose the function", () => {
  const D = disc("1/3 + 1/5i", "1/2");
  const cases = [
    ["exp", (z) => exp(z), (p) => exp(p)],
    ["sin", (z) => sin(z), (p) => sin(p)],
    ["cos", (z) => cos(z), (p) => cos(p)],
    ["sinh", (z) => sinh(z), (p) => sinh(p)],
    ["cosh", (z) => cosh(z), (p) => cosh(p)],
    ["tan", (z) => tan(z), (p) => tan(p)],
    ["tanh", (z) => tanh(z), (p) => tanh(p)],
    ["log(z + 2)", (z) => log(z.add(2)), (p) => log(p.add(2)).rect],
    ["sqrt(z + 2)", (z) => sqrt(z.add(2)), (p) => sqrt(p.add(2)).rect],
    ["cube root of z + 3", (z) => root(z.add(3), 3), (p) => root(p.add(3), 3).rect],
    ["(z + 2)^(1/2 + i)", (z) => pow(z.add(2), "1/2 + i"), (p) => pow(p.add(2), "1/2 + i").rect],
    ["2^z", (z) => pow(2, z), (p) => pow(2, p).rect],
    ["exp(z)/(z − 2)", (z) => exp(z).div(z.sub(2)), (p) => exp(p).div(p.sub(2))],
  ];
  for (const [name, f, ref] of cases) {
    test(name, () => {
      // low orders too: there the bound on the dropped tail is what is being tested
      for (const order of [2, 5, 14]) {
        for (const d of [D, disc("-1/4 + 1/3i", "1/8")]) {
          const m = taylorModel(f, d, { order });
          const range = m.range();
          for (const p of pointsOf(d)) {
            const want = ref(p);
            expect(meets(m.at(p), want)).toBe(true);
            expect(Rect.from(range).intersects(want)).toBe(true);
          }
          if (order === 14) expect(m.error.toNumber()).toBeLessThan(1e-2);
        }
      }
    });
  }
  test("a polynomial is enclosed exactly, and tightly", () => {
    const m = taylorModel("z^5 - 3z^2 + (2 - i)z + 1/7", D);
    for (const p of pointsOf(D)) {
      const exact = p.pow(5).sub(p.pow(2).scale(3)).add(complex("2 - i").mul(p)).add(complex("1/7"));
      expect(m.at(p).contains(exact)).toBe(true);
    }
    expect(m.error.toNumber()).toBeLessThan(1e-24);
  });
  test("text and functions build the same model", () => {
    const a = taylorModel("exp(z) - 3z^2 + pi", D);
    const b = taylorModel((z) => exp(z).sub(z.pow(2).scale(3)).add(CTM.constant(D, rect([3.14159265358979, 3.1415926535898], [0, 0]))), D);
    expect(meets(a.at(D.center), b.at(D.center))).toBe(true);
  });
  test("more order means less error; more bits lets it go on shrinking", () => {
    const e = (order, bits) => taylorModel("exp(z)", disc(0, 1), { order, bits }).error.toNumber();
    expect(e(4)).toBeGreaterThan(e(8));
    expect(e(8)).toBeGreaterThan(e(16));
    expect(e(40, 200)).toBeLessThan(1e-45);
  });
  test("coefficients are the Taylor coefficients", () => {
    const m = taylorModel("exp(z)", disc(0, "1/2"), { order: 10 });
    const c = m.coefficients.map((z) => z.toNumbers()[0]);
    expect(c[0]).toBeCloseTo(1, 12);
    expect(c[3]).toBeCloseTo(1 / 6, 12);
    expect(c[5]).toBeCloseTo(1 / 120, 12);
  });
});

describe("what a model refuses", () => {
  const code = (fn) => {
    try {
      fn();
    } catch (e) {
      expect(e).toBeInstanceOf(ComplexError);
      return e.code;
    }
    return "no error";
  };
  test("a divisor that may vanish", () => {
    expect(code(() => taylorModel("1/z", disc("1/4", "1/2")))).toBe("division-by-zero");
    expect(code(() => taylorModel("1/z", disc("1", "1/2")))).toBe("no error");
  });
  test("a logarithm or root at a possible zero", () => {
    expect(code(() => taylorModel("log(z)", disc("1/4", "1/2")))).toBe("branch-point");
    expect(code(() => taylorModel("sqrt(z)", disc(0, "1/2")))).toBe("branch-point");
  });
  test("an exponential that varies too much for the order", () => {
    expect(code(() => taylorModel("exp(z)", disc(0, 40), { order: 6 }))).toBe("precision");
  });
  test("models on different discs do not mix", () => {
    const a = CTM.z(disc(0, 1));
    const b = CTM.z(disc(1, 1));
    expect(code(() => a.add(b))).toBe("bad-argument");
  });
  test("points outside the disc", () => {
    expect(code(() => taylorModel("z", disc(0, 1)).at(2))).toBe("bad-argument");
  });
  test("a zero-radius disc is not a disc", () => {
    expect(() => taylorModel("z", disc(0, 0))).toThrow(RangeError);
  });
});

describe("derivative and integral", () => {
  test("the derivative model encloses the derivative on the smaller disc", () => {
    const m = taylorModel("exp(z) + z^3", disc("1/4", 1), { order: 16 });
    const d = m.derivative("1/2");
    for (const p of pointsOf(d.disc)) {
      expect(meets(d.at(p), exp(p).add(p.pow(2).scale(3)))).toBe(true);
    }
    expect(() => m.derivative(1)).toThrow(RangeError);
  });
  test("the integral of a polynomial is exact to the working precision", () => {
    const m = taylorModel("3z^2 + 2z", disc(0, 2));
    const v = m.integral("-1", "1 + i");
    // z³ + z² from −1 to 1 + i: (−2 + 2i) + 2i − (−1 + 1) = −2 + 4i
    expect(v.contains(complex("-2 + 4i"))).toBe(true);
    expect(v.size.toNumber()).toBeLessThan(1e-24);
  });
  test("the integral of exp is exp at the ends", () => {
    const m = taylorModel("exp(z)", disc(0, 1), { order: 20 });
    const v = m.integral("-1/2", "1/2 + 1/2i");
    expect(meets(v, exp("1/2 + 1/2i").sub(exp("-1/2")))).toBe(true);
    expect(v.size.toNumber()).toBeLessThan(1e-15);
  });
});

describe("branches and zero counts", () => {
  const code = (fn) => {
    try {
      fn();
    } catch (e) {
      return e.code;
    }
    return "no error";
  };
  test("a principal branch is refused on a disc that reaches its cut", () => {
    expect(code(() => taylorModel("log(z)", disc("-1", "1/4")))).toBe("branch-cut");
    expect(code(() => taylorModel("sqrt(z)", disc("-1 + 1/8i", "1/4")))).toBe("branch-cut");
    expect(code(() => taylorModel("log(z)", disc("-1 + i", "1/4")))).toBe("no error");
    // the same logarithm with its cut turned to the positive axis
    const m = taylorModel("log(-z) + i*pi", disc("-1", "1/4"));
    expect(Math.abs(m.at(-1).center.toNumbers()[1] - Math.PI)).toBeLessThan(1e-9);
  });
  test("a numbered branch of log adds 2πi", () => {
    const a = taylorModel((z) => log(z), disc("2 + i", "1/4"));
    const b = taylorModel((z) => log(z, { branch: 1 }), disc("2 + i", "1/4"));
    const d = b.at("2 + i").sub(a.at("2 + i")).center.toNumbers();
    expect(d[0]).toBeCloseTo(0, 12);
    expect(d[1]).toBeCloseTo(2 * Math.PI, 12);
  });
  test("zeroCount is Rouché's theorem on the coefficients", () => {
    expect(taylorModel("z^3", disc(0, 1)).zeroCount()).toBe(3);
    expect(taylorModel("z^3 + 1/100", disc(0, 1)).zeroCount()).toBe(3);
    expect(taylorModel("z^2 - 2", disc("3/2", "1/4")).zeroCount()).toBe(1);
    expect(taylorModel("exp(z)", disc(0, "1/4")).zeroCount()).toBe(0);
    // zeros at ±1 on the rim: nothing can be said
    expect(taylorModel("z^2 - 1", disc(0, 1)).zeroCount()).toBe(null);
  });
});
