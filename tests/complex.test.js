import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval } from "@ratmath/core";
import * as C from "../index.js";

const { complex, rect, disc } = C;
const s = (x) => String(x);
/** The centre of an enclosure as doubles, and its half-size. */
const mid = (r) => r.center.toNumbers();
const near = (r, [re, im], tol = 1e-12) => {
  const [a, b] = mid(r);
  expect(Math.abs(a - re)).toBeLessThan(tol);
  expect(Math.abs(b - im)).toBeLessThan(tol);
};

describe("complex rationals", () => {
  test("reading and printing", () => {
    expect(s(complex("3 - 4i"))).toBe("3 - 4i");
    expect(s(complex("1/2 + 3/4i"))).toBe("1/2 + 3/4i");
    expect(s(complex("(1 + 2i)/(3 - i)"))).toBe("1/10 + 7/10i");
    expect(s(complex("(1+i)^8"))).toBe("16");
    expect(s(complex("2i"))).toBe("2i");
    expect(s(complex("-i"))).toBe("-i");
    expect(s(complex("0.#3 + 1..1/2 i"))).toBe("1/3 + 3/2i");
    expect(s(complex(3, -4))).toBe("3 - 4i");
    expect(s(complex([1, 2]))).toBe("1 + 2i");
    expect(s(complex(0.1))).toBe("1/10");
    expect(complex("1/2 - i").toTeX()).toBe("\\frac{1}{2} - i");
    expect(() => complex("1 +")).toThrow(SyntaxError);
    expect(() => complex("2^i")).toThrow(SyntaxError);
  });
  test("exact arithmetic", () => {
    const z = complex("3 + 4i");
    expect(s(z.mul(z.conj()))).toBe("25");
    expect(s(z.norm())).toBe("25");
    expect(s(z.exactModulus())).toBe("5");
    expect(complex("1 + i").exactModulus()).toBe(null);
    expect(s(z.inv())).toBe("3/25 - 4/25i");
    expect(s(z.div("1 - 2i"))).toBe("-1 + 2i");
    expect(s(z.pow(-2))).toBe("-7/625 - 24/625i");
    expect(s(complex("i").pow(1001))).toBe("i");
    expect(z.equals("3 + 4i")).toBe(true);
    expect(complex("2i").isImaginary()).toBe(true);
    expect(z.isGaussianInteger()).toBe(true);
    expect(() => z.div(0)).toThrow(C.ComplexError);
  });
  test("exact square roots", () => {
    expect(complex("3 + 4i").exactSqrt().map(s)).toEqual(["2 + i", "-2 - i"]);
    expect(complex("-4").exactSqrt().map(s)).toEqual(["2i", "-2i"]);
    expect(complex("2").exactSqrt()).toBe(null);
    expect(s(complex("-5 - 12i").exactSqrt()[0])).toBe("2 - 3i");
  });
  test("the field object", () => {
    const { CQ } = C;
    expect(s(CQ.mul(CQ.from("1 + i"), CQ.from("1 - i")))).toBe("2");
    expect(CQ.sign(CQ.zero)).toBe(0);
    expect(CQ.equals(CQ.inv("i"), "-i")).toBe(true);
  });
});

describe("rectangles and discs", () => {
  test("rectangle arithmetic contains the pointwise results", () => {
    const a = rect([1, 2], [0, 1]);
    const b = rect(["-1/2", "1/2"], [2, 3]);
    const pts = (r) => [...r.corners(), r.center];
    for (const [name, f] of [["add", (x, y) => x.add(y)], ["sub", (x, y) => x.sub(y)], ["mul", (x, y) => x.mul(y)], ["div", (x, y) => x.div(y)]]) {
      const out = a[name](b);
      for (const p of pts(a)) for (const q of pts(b)) expect(out.contains(f(p, q))).toBe(true);
    }
    expect(a.pow(3).contains(complex("3/2 + 1/2i").pow(3))).toBe(true);
    expect(() => a.div(rect([-1, 1], [-1, 1]))).toThrow(C.ComplexError);
  });
  test("shape", () => {
    const r = rect("1 + i", "1/10");
    expect(s(r.center)).toBe("1 + i");
    expect(s(r.width)).toBe("1/5");
    expect(r.contains("21/20 + 19/20i")).toBe(true);
    expect(r.contains("2")).toBe(false);
    expect(r.isPoint).toBe(false);
    expect(rect("2i").isPoint).toBe(true);
    expect(r.split().length).toBe(4);
    expect(r.split().every((q) => r.contains(q))).toBe(true);
    expect(r.intersect(rect("3", "1/10"))).toBe(null);
    expect(r.hull(rect("3", "1/10")).contains("2 + 1/2i")).toBe(true);
    expect(r.real).toBeInstanceOf(RationalInterval);
    expect(r.toDisc().contains(r)).toBe(true);
    expect(s(rect("3 - 4i"))).toBe("3 - 4i");
  });
  test("disc arithmetic", () => {
    const a = disc("1 + i", "1/10");
    const b = disc("2", "1/10");
    const p = a.mul(b);
    expect(s(p.center)).toBe("2 + 2i");
    expect(p.contains(complex("11/10 + i").mul("21/10"))).toBe(true);
    // 1/z of a disc is a disc, exactly
    const inv = disc("2", 1).inv();
    expect(s(inv.center)).toBe("2/3");
    expect(s(inv.radius)).toBe("1/3");
    expect(disc(0, 1).contains("3/5 + 4/5i")).toBe(true);
    expect(disc(0, 1).contains("3/5 + 5/6i")).toBe(false);
    expect(disc(0, 1).intersects(disc(2, 1))).toBe(true);
    expect(disc(0, 1).intersects(disc("2 + 1/100i", 1))).toBe(false);
    const h = disc(1, "1/2").hull(disc(3, "1/4"));
    expect(h.contains(disc(1, "1/2")) && h.contains(disc(3, "1/4"))).toBe(true);
    expect(() => disc(0, 1).inv()).toThrow(C.ComplexError);
  });
});

describe("polar form", () => {
  test("modulus and argument", () => {
    const p = C.polar("3 + 4i");
    expect(s(p.exactModulus)).toBe("5");
    expect(s(p.normSquared)).toBe("25");
    expect(p.exactTurn).toBe(null);
    expect(p.modulus.low.equals(new Rational(5n))).toBe(true);
    const q = C.polar("1 + i");
    expect(q.exactModulus).toBe(null);
    expect(s(q.exactTurn)).toBe("1/8");
    expect(q.modulus.low.toNumber()).toBeCloseTo(Math.SQRT2, 14);
    expect(q.argument.low.toNumber()).toBeCloseTo(Math.PI / 4, 14);
    expect(q.argument.high.subtract(q.argument.low).toNumber()).toBeLessThan(1e-18);
    expect(s(C.exactTurn("-i"))).toBe("-1/4");
    expect(s(C.exactTurn("-3"))).toBe("1/2");
    expect(C.argument("-1 - i").low.toNumber()).toBeCloseTo((-3 * Math.PI) / 4, 14);
  });
  test("back from polar form", () => {
    expect(s(C.fromPolar(2, { turns: "1/4" }))).toBe("2i");
    expect(s(C.fromPolar(3, { turns: "1/2" }))).toBe("-3");
    near(C.fromPolar(2, { turns: "1/8" }), [Math.SQRT2, Math.SQRT2]);
    near(C.fromPolar(1, "1/2"), [Math.cos(0.5), Math.sin(0.5)]);
    expect(C.fromPolar(5, C.argument("3 + 4i")).contains("3 + 4i")).toBe(true);
  });
});

describe("functions", () => {
  const zs = ["1 + i", "-2 + 1/3i", "1/4 - 3i", "-1/2 - 1/2i", "3", "2i"];
  const cx = (z) => complex(z).toNumbers();
  test("exp, sin, cos and the hyperbolic ones enclose the true value", () => {
    for (const z of zs) {
      const [x, y] = cx(z);
      near(C.exp(z), [Math.exp(x) * Math.cos(y), Math.exp(x) * Math.sin(y)], 1e-10);
      near(C.sin(z), [Math.sin(x) * Math.cosh(y), Math.cos(x) * Math.sinh(y)], 1e-10);
      near(C.cos(z), [Math.cos(x) * Math.cosh(y), -Math.sin(x) * Math.sinh(y)], 1e-10);
      near(C.sinh(z), [Math.sinh(x) * Math.cos(y), Math.cosh(x) * Math.sin(y)], 1e-10);
      near(C.cosh(z), [Math.cosh(x) * Math.cos(y), Math.sinh(x) * Math.sin(y)], 1e-10);
    }
  });
  test("identities hold as enclosures", () => {
    for (const z of zs) {
      const one = C.sin(z).pow(2).add(C.cos(z).pow(2));
      expect(one.contains(1)).toBe(true);
      expect(one.size.toNumber()).toBeLessThan(1e-12);
      expect(C.tan(z).mul(C.cos(z)).intersects(C.sin(z))).toBe(true);
      expect(C.tanh(z).mul(C.cosh(z)).intersects(C.sinh(z))).toBe(true);
      expect(C.exp(C.log(z, { side: "above" }).rect).contains(z)).toBe(true);
    }
  });
  test("precision is a choice", () => {
    const r = C.exp("1 + i", 300);
    expect(r.size.toNumber()).toBeLessThan(1e-85);
    expect(C.exp("1 + i", 20).contains(r)).toBe(true);
    expect(C.exp("1 + i", { bits: 20 }).size.toNumber()).toBeLessThan(1e-5);
  });
  test("a wide input gives a wide answer that still contains every value", () => {
    const box = rect("1 + i", "1/10");
    const out = C.exp(box);
    for (const c of [...box.corners(), box.center]) expect(out.contains(C.exp(c))).toBe(true);
    expect(() => C.tan(rect([1, 2], [-1, 1]))).toThrow(C.ComplexError); // contains π/2
  });
});

describe("branches", () => {
  test("log away from its cut", () => {
    const l = C.log("1 + i");
    expect(l.straddles).toBe(false);
    expect(l.branch).toBe("principal");
    near(l.value, [Math.log(2) / 2, Math.PI / 4]);
    near(C.log("i", { branch: 1 }).value, [0, Math.PI / 2 + 2 * Math.PI]);
    expect(C.log("i", { branch: 1 }).branch).toBe(1);
  });
  test("log on its cut reports both sides and the convention", () => {
    const l = C.log("-1");
    expect(l.straddles).toBe(true);
    expect(l.value).toBe(null);
    expect(l.candidates.map((c) => c.side)).toEqual(["above", "below"]);
    near(l.candidates[0].value, [0, Math.PI]);
    near(l.candidates[1].value, [0, -Math.PI]);
    expect(l.convention.side).toBe("above");
    expect(() => l.rect).toThrow(/branch cut/);
    expect(l.hull.contains(l.candidates[0].value) && l.hull.contains(l.candidates[1].value)).toBe(true);
    near(C.log("-1", { side: "below" }).rect, [0, -Math.PI]);
    expect(String(l)).toMatch(/^on the cut/);
  });
  test("an enclosure across the cut has no convention", () => {
    const l = C.log(rect("-2", "1/10"));
    expect(l.straddles).toBe(true);
    expect(l.convention).toBe(null);
    expect(l.candidates.length).toBe(2);
    expect(C.log(rect("-2 + i", "1/10")).straddles).toBe(false);
    expect(() => C.log(rect("0", "1/10"), { side: "above" })).toThrow(C.ComplexError);
  });
  test("square roots", () => {
    expect(s(C.sqrt("3 + 4i").value)).toBe("2 + i");
    expect(s(C.sqrt("3 + 4i", { branch: 1 }).value)).toBe("-2 - i");
    const r = C.sqrt("-4");
    expect(r.straddles).toBe(true);
    expect(r.candidates.map((c) => s(c.value))).toEqual(["2i", "-2i"]);
    expect(s(C.sqrt("-4", { side: "below" }).value)).toBe("-2i");
    near(C.sqrt("2").value, [Math.SQRT2, 0]);
    near(C.sqrt("i").value, [Math.SQRT1_2, Math.SQRT1_2]);
    expect(C.sqrt("-2").straddles).toBe(true);
  });
  test("all the roots", () => {
    const rs = C.roots("8", 3);
    expect(rs.length).toBe(3);
    near(rs[0], [2, 0]);
    near(rs[1], [-1, Math.sqrt(3)]);
    near(rs[2], [-1, -Math.sqrt(3)]);
    for (const r of rs) expect(r.pow(3).contains(8)).toBe(true);
    const u = C.rootsOfUnity(8);
    expect([0, 2, 4, 6].map((k) => s(u[k]))).toEqual(["1", "i", "-1", "-i"]);
    near(u[1], [Math.SQRT1_2, Math.SQRT1_2]);
    expect(() => C.root("2", 0)).toThrow(RangeError);
  });
  test("powers", () => {
    near(C.pow("i", "i").value, [Math.exp(-Math.PI / 2), 0]);
    expect(s(C.pow("1 + i", 4).value)).toBe("-4");
    near(C.pow("2", "1/2").value, [Math.SQRT2, 0]);
    expect(s(C.pow("4", "3/2").value)).toBe("8");
    expect(C.pow("-1", "1/2").straddles).toBe(true);
    near(C.pow("i", "i", { branch: 1 }).value, [Math.exp(-Math.PI / 2 - 2 * Math.PI), 0]);
  });
  test("inverse functions and their cuts", () => {
    near(C.asin("1/2").value, [Math.PI / 6, 0]);
    near(C.acos("1/2").value, [Math.PI / 3, 0]);
    near(C.atan("1").value, [Math.PI / 4, 0]);
    near(C.asinh("1").value, [Math.asinh(1), 0]);
    near(C.acosh("2").value, [Math.acosh(2), 0]);
    near(C.atanh("1/2").value, [Math.atanh(0.5), 0]);
    // on the cuts
    const a = C.asin("2");
    expect(a.straddles).toBe(true);
    const im = Math.log(2 + Math.sqrt(3));
    near(a.candidates.find((c) => c.side === "above").value, [Math.PI / 2, im]);
    near(a.candidates.find((c) => c.side === "below").value, [Math.PI / 2, -im]);
    expect(a.convention.side).toBe("below");
    expect(C.asin("-2").convention.side).toBe("above");
    const t = C.atan("2i");
    expect(t.straddles).toBe(true);
    expect(t.candidates.map((c) => c.side)).toEqual(["right", "left"]);
    expect(C.acosh("-3").straddles).toBe(true);
    near(C.acosh("-3", { side: "above" }).rect, [Math.acosh(3), Math.PI]);
    expect(C.atanh("2").straddles).toBe(true);
    expect(C.asinh("2i").straddles).toBe(true);
    expect(C.acos("3").straddles).toBe(true);
    // each inverse undoes its function off the cut
    for (const z of ["1 + i", "-1/3 + 2i", "1/2 - 1/4i"]) {
      expect(C.sin(C.asin(z).rect).contains(z)).toBe(true);
      expect(C.cos(C.acos(z).rect).contains(z)).toBe(true);
      expect(C.tan(C.atan(z).rect).contains(z)).toBe(true);
      expect(C.sinh(C.asinh(z).rect).contains(z)).toBe(true);
      expect(C.cosh(C.acosh(z).rect).contains(z)).toBe(true);
      expect(C.tanh(C.atanh(z).rect).contains(z)).toBe(true);
    }
    expect(() => C.asin("2", { branch: 1 })).toThrow(C.ComplexError);
  });
  test("the argument as a branched result", () => {
    const a = C.arg("-1");
    expect(a.straddles).toBe(true);
    near(a.candidates[0].value, [Math.PI, 0]);
    near(C.arg("i").value, [Math.PI / 2, 0]);
  });
});

describe("continuation", () => {
  test("winding numbers are exact", () => {
    expect(C.winding(["1", "i", "-1", "-i"])).toBe(1);
    expect(C.winding(["1", "-i", "-1", "i"])).toBe(-1);
    expect(C.winding(["1", "2", "2 + i"])).toBe(0);
    expect(C.winding(["2", "2i", "-2", "-2i"], "1/2")).toBe(1);
    expect(C.winding(["1", "i", "-1", "-i", "1", "i", "-1", "-i"])).toBe(2);
    expect(() => C.winding(["1", "-1"])).toThrow(C.ComplexError);
  });
  test("log carried once round the origin gains 2πi", () => {
    const r = C.continueAlong("log", ["1", "i", "-1", "-i", "1"]);
    expect(r.turns).toBe(1);
    near(r.value, [0, 2 * Math.PI]);
    expect(r.values.length).toBe(5);
    near(r.values[2], [0, Math.PI]);
  });
  test("a square root comes back with the other sign, a cube root needs three turns", () => {
    const loop = ["1", "i", "-1", "-i", "1"];
    near(C.continueAlong({ root: 2 }, loop).value, [-1, 0]);
    const three = [...loop, ...loop.slice(1), ...loop.slice(1)];
    near(C.continueAlong({ root: 3 }, loop).value, [-0.5, Math.sqrt(3) / 2]);
    near(C.continueAlong({ root: 3 }, three).value, [1, 0]);
    near(C.continueAlong({ pow: "1/2" }, loop).value, [-1, 0]);
    near(C.continueAlong("log", ["2", "1 + i", "0", "1 - i", "2"], { about: 1 }).value, [0, 2 * Math.PI]);
    expect(() => C.continueAlong("log", ["1", "-1"])).toThrow(C.ComplexError);
  });
});

test("refined gives a quick value, then the certified one", () => {
  const r = C.refined("exp", "1 + i", { bits: 100 });
  expect(r.quick.value[0]).toBeCloseTo(Math.E * Math.cos(1), 6);
  expect(r.quick.trust).toBe("good");
  const exact = r.exact();
  expect(exact.size.toNumber()).toBeLessThan(1e-28);
  expect(C.refined("log", "-1").quick.trust).toBe("unsure");
  expect(() => C.refined("nope", 1)).toThrow(RangeError);
});
