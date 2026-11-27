import { expect, test, describe } from "bun:test";
import { Contour, contourIntegral, residueSum, argumentPrinciple, countZeros, zeros, complex, rect, disc, exp, sqrt, Rect, ComplexError } from "../index.js";

const near = (r, re, im, eps = 1e-9) => {
  const [x, y] = r.center.toNumbers();
  expect(Math.abs(x - re)).toBeLessThan(eps);
  expect(Math.abs(y - im)).toBeLessThan(eps);
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

describe("contours", () => {
  test("polygons, rectangles and circles are closed; paths are not", () => {
    expect(Contour.polygon(["0", "1", "i"]).closed).toBe(true);
    expect(Contour.path(["0", "1", "i"]).closed).toBe(false);
    expect(Contour.rectangle(rect([0, 1], [0, 2])).pieces.length).toBe(4);
    const c = Contour.circle("1 + i", 2);
    expect(c.closed).toBe(true);
    expect(c.pieces.length).toBe(4);
    // every joint is an exact point of the circle
    for (const p of c.pieces) expect(p.a.sub(complex("1 + i")).norm().toString()).toBe("4");
  });
  test("arcs go the way asked", () => {
    const short = Contour.arc(0, 1, "i");
    const long = Contour.arc(0, 1, "i", { clockwise: true });
    expect(short.pieces.length).toBe(1);
    expect(long.pieces.length).toBeGreaterThanOrEqual(3);
    expect(short.join(long.reversed()).closed).toBe(true);
    expect(() => Contour.arc(0, 1, 2)).toThrow(RangeError);
  });
  test("a broken contour is refused", () => {
    expect(() => Contour.path(["0", "1"]).join(Contour.path(["2", "3"]))).toThrow(RangeError);
  });
  test("winding numbers", () => {
    expect(Contour.circle("1/3", 2).winding(0)).toBe(1);
    expect(Contour.circle("1/3", 2, { clockwise: true }).winding(0)).toBe(-1);
    expect(Contour.circle("1/3", 2).winding(5)).toBe(0);
    const twice = Contour.circle(0, 1).join(Contour.circle(0, 1));
    expect(twice.winding("1/2")).toBe(2);
    expect(code(() => Contour.circle(0, 1).winding(1))).toBe("zero-on-contour");
  });
  test("sample gives points for drawing", () => {
    const pts = Contour.circle(0, 1).sample(4);
    expect(pts.length).toBe(17);
    for (const [x, y] of pts) expect(Math.hypot(x, y)).toBeCloseTo(1, 10);
  });
});

describe("contour integrals", () => {
  test("a polynomial along a path", () => {
    const v = contourIntegral("z^2", Contour.path(["0", "1", "1 + i"]));
    expect(v.contains(complex("-2/3 + 2/3i"))).toBe(true);
    expect(v.size.toNumber()).toBeLessThan(1e-20);
  });
  test("an analytic function integrates to the same thing along different paths", () => {
    const a = contourIntegral("exp(z) * z", Contour.path(["0", "2", "2 + i"]));
    const b = contourIntegral("exp(z) * z", Contour.path(["0", "i", "2 + i"]));
    const c = contourIntegral("exp(z) * z", Contour.arc("1 + 1/2i", "0", "2 + i"));
    expect(a.intersects(b)).toBe(true);
    expect(a.intersects(c)).toBe(true);
    // (z − 1)e^z from 0 to 2 + i
    const w = complex("2 + i");
    expect(a.intersects(exp(w).mul(w.sub(1)).add(1))).toBe(true);
  });
  test("1/z round the origin is 2πi, and −2πi the other way", () => {
    const v = contourIntegral("1/z", Contour.circle(0, 1));
    near(v, 0, 2 * Math.PI);
    expect(v.size.toNumber()).toBeLessThan(1e-11);
    near(contourIntegral("1/z", Contour.circle(0, 1, { clockwise: true })), 0, -2 * Math.PI);
    near(contourIntegral("1/z", rect([-1, 2], [-1, 1])), 0, 2 * Math.PI);
    near(contourIntegral("1/z", Contour.circle(3, 1)), 0, 0);
  });
  test("a tighter tolerance is met", () => {
    const v = contourIntegral("1/z", Contour.circle(0, 1), { tolerance: "1/1000000000000000000000000", order: 20, bits: 128 });
    expect(v.size.toNumber()).toBeLessThan(1e-24);
    near(v, 0, 2 * Math.PI, 1e-14);
  });
  test("Cauchy's integral formula", () => {
    // ∮ e^z/(z − 1/2) dz = 2πi·e^(1/2)
    const v = contourIntegral("exp(z)/(z - 1/2)", Contour.circle(0, 1));
    near(v, 0, 2 * Math.PI * Math.exp(0.5));
  });
  test("residues", () => {
    near(residueSum("exp(z)/z^3", Contour.circle(0, 1)), 0.5, 0);
    near(residueSum("1/(z^2 + 1)", rect([-1, 1], [0, 2])), 0, -0.5);
    near(residueSum("1/(z^2 + 1)", Contour.circle(0, 2)), 0, 0);
    expect(code(() => residueSum("1/z", Contour.path(["1", "2"])))).toBe("bad-argument");
  });
  test("a function given as a JavaScript function, with a branch", () => {
    // sqrt(z + 3) is analytic on the unit disc
    const v = contourIntegral((z) => sqrt(z.add(3)), Contour.circle(0, 1));
    near(v, 0, 0);
  });
  test("a pole on the contour is refused", () => {
    expect(code(() => contourIntegral("1/(z - 1)", Contour.circle(0, 1), { maxDepth: 12 }))).toBe("singularity");
  });
});

describe("the argument principle", () => {
  test("counts zeros", () => {
    expect(argumentPrinciple("z^3 - 1", Contour.circle(0, 2))).toBe(3);
    expect(argumentPrinciple("z^3 - 1", Contour.circle(1, "1/2"))).toBe(1);
    expect(countZeros("exp(z) - 3z", rect([-1, 3], [-2, 2]))).toBe(2);
    expect(countZeros("z^5", disc(0, 1))).toBe(5);
  });
  test("counts zeros minus poles", () => {
    // tan has a zero at 0 and poles at ±π/2 inside |z| = 2
    expect(argumentPrinciple("tan(z)", Contour.circle(0, 2))).toBe(-1);
    expect(argumentPrinciple("(z - 1)^2/(z + 1)^3", Contour.circle(0, 3))).toBe(-1);
  });
  test("a clockwise contour counts negatively", () => {
    expect(argumentPrinciple("z^3 - 1", Contour.circle(0, 2, { clockwise: true }))).toBe(-3);
  });
  test("a zero on the contour is refused", () => {
    expect(code(() => argumentPrinciple("z - 1", Contour.circle(0, 1)))).toBe("zero-on-contour");
    expect(code(() => argumentPrinciple("z", Contour.path(["1", "2"])))).toBe("bad-argument");
  });
});

describe("zeros", () => {
  test("simple zeros, each proved and refined", () => {
    const r = zeros("z^2 - 2", rect([-2, 2], [-1, 1]));
    expect(r.complete).toBe(true);
    expect(r.count).toBe(2);
    for (const z of r.zeros) {
      expect(z.simple).toBe(true);
      expect(z.multiplicity).toBe(1);
      expect(z.enclosure.size.toNumber()).toBeLessThan(1e-18);
      // √2 is really inside: x² − 2 changes sign across the enclosure's real part
      const { lo, hi } = z.enclosure.re;
      const [a, b] = [lo.multiply(lo).toNumber() - 2, hi.multiply(hi).toNumber() - 2];
      expect(lo.multiply(lo).compareTo(hi.multiply(hi)) !== 0).toBe(true);
      expect(Math.abs(a)).toBeLessThan(1e-17);
      expect(Math.abs(b)).toBeLessThan(1e-17);
    }
  });
  test("exact rational zeros are inside their enclosures", () => {
    const r = zeros("(z - 1)*(z - 2)*(z - 3)*(z - 1/2 - 1/3i)", rect([0, 4], [-1, 1]));
    expect(r.count).toBe(4);
    for (const w of ["1", "2", "3", "1/2 + 1/3i"]) expect(r.zeros.some((z) => z.enclosure.contains(complex(w)))).toBe(true);
  });
  test("a transcendental equation, checked against the count on the boundary", () => {
    const R = rect([-10, 10], [-10, 10]);
    const r = zeros("cos(z) - z", R);
    expect(r.complete).toBe(true);
    expect(r.count).toBe(7);
    expect(countZeros("cos(z) - z", R)).toBe(7);
    // the real one is the Dottie number
    expect(r.zeros.some((z) => Math.abs(z.enclosure.center.toNumbers()[0] - 0.7390851332151607) < 1e-15 && z.enclosure.im.lo.numerator <= 0n && z.enclosure.im.hi.numerator >= 0n)).toBe(true);
  });
  test("the zeros of exp(z) − 1 are the multiples of 2πi", () => {
    const r = zeros("exp(z) - 1", rect([-1, 1], [-10, 10]));
    expect(r.count).toBe(3);
    expect(r.zeros.map((z) => Math.round(z.enclosure.center.toNumbers()[1] / (2 * Math.PI))).sort((a, b) => a - b)).toEqual([-1, 0, 1]);
  });
  test("a double zero is reported with its multiplicity", () => {
    const r = zeros("(z - 1/3)^2 * (z + 1)", rect([-2, 2], [-1, 1]));
    expect(r.complete).toBe(true);
    expect(r.count).toBe(3);
    const double = r.zeros.find((z) => !z.simple);
    expect(double.multiplicity).toBe(2);
    expect(double.enclosure.contains(complex("1/3"))).toBe(true);
    expect(double.enclosure.size.toNumber()).toBeLessThan(1e-5);
    expect(r.zeros.find((z) => z.simple).enclosure.contains(complex(-1))).toBe(true);
  });
  test("a triple zero at a corner of the cells", () => {
    const r = zeros("z^3", rect([-1, 1], [-1, 1]));
    expect(r.zeros.length).toBe(1);
    expect(r.zeros[0].multiplicity).toBe(3);
  });
  test("two zeros closer than the tolerance are counted together, and apart at a finer one", () => {
    const f = "(z - 1/3)*(z - 1/3 - 1/1000000000)";
    const coarse = zeros(f, rect([-1, 1], [-1, 1]));
    expect(coarse.count).toBe(2);
    expect(coarse.zeros.length).toBe(1);
    expect(coarse.zeros[0].simple).toBe(false);
    const fine = zeros(f, rect([-1, 1], [-1, 1]), { tolerance: "1/1000000000000" });
    expect(fine.zeros.length).toBe(2);
    expect(fine.zeros.every((z) => z.simple)).toBe(true);
  });
  test("a zero on the edge is flagged, not counted", () => {
    const r = zeros("z*(z - 1/2)", rect([0, 1], [-1, 1]));
    expect(r.count).toBe(1);
    expect(r.zeros.length).toBe(2);
    expect(r.zeros.find((z) => z.boundary).enclosure.contains(complex(0))).toBe(true);
  });
  test("a pole in the region leaves it unresolved, and says where", () => {
    const r = zeros("tan(z)", rect([-1, 2], [-1, 1]), { tolerance: "1/1000" });
    expect(r.complete).toBe(false);
    expect(r.zeros.length).toBe(1);
    expect(r.unresolved.length).toBeGreaterThan(0);
    for (const u of r.unresolved) expect(Math.abs(u.region.center.toNumbers()[0] - Math.PI / 2)).toBeLessThan(0.01);
  });
  test("the zero function runs out of budget instead of running forever", () => {
    const r = zeros("z - z", rect([0, 1], [-1, 1]), { maxCells: 300 });
    expect(r.complete).toBe(false);
    expect(r.zeros.length).toBe(0);
  });
  test("more bits, more digits", () => {
    const r = zeros("z^2 - 2", rect([1, 2], [-1, 1]), { bits: 200 });
    expect(r.zeros[0].enclosure.size.toNumber()).toBeLessThan(1e-58);
  });
  test("a JavaScript function works too", () => {
    const r = zeros((z) => z.pow(5).sub(z).add(1), rect([-2, 2], [-2, 2]));
    expect(r.count).toBe(5);
  });
});
