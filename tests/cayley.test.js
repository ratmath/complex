import { describe, expect, test } from "bun:test";
import * as C from "../cayley.js";
import { CQ, complex } from "../index.js";

const { quaternions: H, octonions: O, sedenions: S, quaternion, octonion } = C;
const s = (x) => String(x);

describe("the doubling construction", () => {
  test("Hamilton's rules", () => {
    const [one, i, j, k] = H.units;
    expect(s(i.mul(j))).toBe("k");
    expect(s(j.mul(k))).toBe("i");
    expect(s(k.mul(i))).toBe("j");
    expect(s(j.mul(i))).toBe("-k");
    for (const u of [i, j, k]) expect(s(u.mul(u))).toBe("-1");
    expect(s(i.mul(j).mul(k))).toBe("-1");
    expect(one.equals(1)).toBe(true);
  });
  test("each level loses a law", () => {
    const p = [0, 1, 2, 3, 4, 5].map((l) => C.cayleyDickson(l).properties());
    expect(p.map((x) => x.dimension)).toEqual([1, 2, 4, 8, 16, 32]);
    expect(p.map((x) => x.commutative)).toEqual([true, true, false, false, false, false]);
    expect(p.map((x) => x.associative)).toEqual([true, true, true, false, false, false]);
    expect(p.map((x) => x.alternative)).toEqual([true, true, true, true, false, false]);
    expect(p.map((x) => x.flexible)).toEqual([true, true, true, true, true, true]);
    expect(p.map((x) => x.zeroDivisors === null)).toEqual([true, true, true, true, false, false]);
    const [a, b] = p[4].zeroDivisors;
    expect(a.isZero() || b.isZero()).toBe(false);
    expect(a.mul(b).isZero()).toBe(true);
    expect(s(H)).toBe("quaternions over Q");
    expect(() => C.cayleyDickson(9)).toThrow(RangeError);
    expect(C.cayleyDickson(8).dim).toBe(256);
  });
  test("level 1 is the complex numbers", () => {
    const Cx = C.cayleyDickson(1);
    const z = Cx.of(3, 4).mul(Cx.of(1, -2));
    const w = complex("3 + 4i").mul("1 - 2i");
    expect(z.coeffs.map(s)).toEqual([s(w.re), s(w.im)]);
    expect(s(Cx.of(3, 4))).toBe("3 + 4i");
  });
  test("reading and printing", () => {
    expect(s(quaternion("1 + 2i - 3/4 k"))).toBe("1 + 2i - 3/4k");
    expect(s(quaternion(1, 2, 3, 4))).toBe("1 + 2i + 3j + 4k");
    expect(s(quaternion("-j"))).toBe("-j");
    expect(s(octonion("e1 + 2 e7"))).toBe("e1 + 2 e7");
    expect(s(octonion(1, 0, 0, -1))).toBe("1 - e3");
    expect(s(H.zero)).toBe("0");
    expect(() => H.parse("1 + e5")).toThrow(SyntaxError);
    expect(() => H.parse("x")).toThrow(SyntaxError);
    expect(() => H.of(1, 2, 3, 4, 5)).toThrow(RangeError);
    expect(s(O.of(quaternion(1, 2, 3, 4)))).toBe("1 + 2 e1 + 3 e2 + 4 e3");
  });
  test("norm, conjugate, inverse", () => {
    const q = quaternion(1, 2, 3, 4);
    expect(s(q.norm())).toBe("30");
    expect(s(q.conj())).toBe("1 - 2i - 3j - 4k");
    expect(s(q.mul(q.conj()))).toBe("30");
    expect(s(q.mul(q.inv()))).toBe("1");
    expect(s(q.real)).toBe("1");
    expect(s(q.imaginary)).toBe("2i + 3j + 4k");
    expect(s(q.scale("1/2"))).toBe("1/2 + i + 3/2j + 2k");
    expect(q.pow(3).equals(q.mul(q).mul(q))).toBe(true);
    expect(q.pow(-1).equals(q.inv())).toBe(true);
    expect(() => H.zero.inv()).toThrow(/inverse/);
    const p = quaternion(2, -1, 0, 5);
    expect(q.mul(p).norm().equals(q.norm().multiply(p.norm()))).toBe(true);
    expect(q.commutator(p).isZero()).toBe(false);
    expect(q.associator(p, quaternion(0, 1, 1, 0)).isZero()).toBe(true);
    expect(q.rightDivide(p).mul(p).equals(q)).toBe(true);
    expect(p.mul(q.leftDivide(p)).equals(q)).toBe(true);
  });
  test("octonions: a norm that multiplies, without associativity", () => {
    const a = O.of(1, 2, 3, 4, 5, 6, 7, 8);
    const b = O.of(2, -1, 0, 3, 1, 1, -2, 5);
    expect(a.mul(b).norm().equals(a.norm().multiply(b.norm()))).toBe(true);
    expect(O.basis(1).associator(O.basis(2), O.basis(4)).isZero()).toBe(false);
    // alternative: (aa)b = a(ab)
    expect(a.mul(a).mul(b).equals(a.mul(a.mul(b)))).toBe(true);
    expect(a.rightDivide(b).mul(b).equals(a)).toBe(true);
    expect(b.mul(a.leftDivide(b)).equals(a)).toBe(true);
  });
  test("sedenions: inverses but no norm law and no division", () => {
    const a = S.of(...Array.from({ length: 16 }, (_, n) => n + 1));
    const b = S.of(...Array.from({ length: 16 }, (_, n) => ((n * 7) % 5) - 2));
    expect(s(a.mul(a.inv()))).toBe("1");
    expect(a.mul(b).norm().equals(a.norm().multiply(b.norm()))).toBe(false);
    expect(a.pow(3).equals(a.mul(a.mul(a)))).toBe(true);
    expect(() => a.rightDivide(b)).toThrow(/not alternative/);
  });
  test("over another field", () => {
    const HC = C.cayleyDickson(2, CQ);
    const q = HC.of("1 + i", "2", "-i", "0");
    expect(q.mul(HC.basis(1)).coeffs.map(s)).toEqual(["-2", "1 + i", "0", "i"]);
    expect(s(q.norm())).toBe("3 + 2i");
    expect(s(q)).toBe("(1 + i) + 2i + (-i)j");
  });
});

describe("rotations", () => {
  const rows = (M) => M.map((r) => r.map(s).join(" "));
  test("a quarter turn about z", () => {
    const q = C.fromAxis([0, 0, 1], 1);
    expect(s(q)).toBe("1 + k");
    expect(C.rotate(q, [1, 0, 0]).map(s)).toEqual(["0", "1", "0"]);
    expect(rows(C.rotationMatrix(q))).toEqual(["0 -1 0", "1 0 0", "0 0 1"]);
    expect(s(C.rotationOf(q).cos)).toBe("0");
    expect(s(C.rotationOf(q).tanHalfSquared)).toBe("1");
  });
  test("every rational quaternion is an exact rotation", () => {
    const q = quaternion(1, 2, 3, 4);
    const M = C.rotationMatrix(q);
    expect(rows(M)).toEqual(["-2/3 2/15 11/15", "2/3 -1/3 2/3", "1/3 14/15 2/15"]);
    // orthogonal with determinant 1
    const dot = (a, b) => a.reduce((t, x, i) => t.add(x.multiply(b[i])), a[0].subtract(a[0]));
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) expect(s(dot(M[i], M[j]))).toBe(i === j ? "1" : "0");
    expect(s(C.fromRotationMatrix(M))).toBe("1 + 2i + 3j + 4k");
    expect(s(C.fromRotationMatrix(C.rotationMatrix("i")))).toBe("i");
    expect(s(C.primitive(quaternion("-1/2", 0, "3/4", 0)))).toBe("2 - 3j");
    expect(() => C.fromRotationMatrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]])).toThrow(/reflection/);
    expect(() => C.fromRotationMatrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])).toThrow(/orthogonal/);
    expect(() => C.rotationMatrix(0)).toThrow(/zero/);
  });
  test("composition is multiplication, in order", () => {
    const a = quaternion(1, 1, 0, 0);
    const b = quaternion(1, 0, 1, 0);
    const v = [1, 2, 3];
    expect(C.rotate(a.mul(b), v).map(s)).toEqual(C.rotate(a, C.rotate(b, v)).map(s));
    expect(C.rotate(a.mul(b), v).map(s)).not.toEqual(C.rotate(b.mul(a), v).map(s));
  });
  test("interpolation stays on rotations", () => {
    const q = C.interpolate(1, "1 + k", "1/2");
    expect(s(q)).toBe("1 + 1/2k");
    expect(C.rotate(q, [1, 0, 0]).map(s)).toEqual(["3/5", "4/5", "0"]);
    // the shorter way round
    expect(s(C.interpolate(1, quaternion(-1, 0, 0, -1), "1/2"))).toBe("1 + 1/2k");
  });
});

describe("integer quaternions", () => {
  const { hurwitz } = C;
  test("Hurwitz quaternions", () => {
    const a = hurwitz(3, 1, -2, 5);
    const b = hurwitz("1/2", "1/2", "1/2", "-1/2");
    expect(a.isLipschitz).toBe(true);
    expect(b.isLipschitz).toBe(false);
    expect(b.norm()).toBe(1n);
    expect(a.norm()).toBe(39n);
    expect(a.mul(b).toQuaternion().equals(a.toQuaternion().mul(b.toQuaternion()))).toBe(true);
    expect(s(a.mul(b))).toBe("9/2 + 1/2i + 7/2j + 5/2k");
    expect(a.toIntegers()).toEqual([3n, 1n, -2n, 5n]);
    expect(() => b.toIntegers()).toThrow(RangeError);
    expect(() => hurwitz("1/2", 1, 0, 0)).toThrow(RangeError);
    expect(() => hurwitz("1/3", 0, 0, 0)).toThrow(RangeError);
    expect(s(a.add(b).sub(b))).toBe(s(a));
    expect(s(a.conj().neg())).toBe("-3 + i - 2j + 5k");
    expect(hurwitz(quaternion(1, 2, 3, 4)).norm()).toBe(30n);
  });
  test("the 24 units", () => {
    const u = C.hurwitzUnits();
    expect(u.length).toBe(24);
    expect(new Set(u.map(s)).size).toBe(24);
    expect(u.every((x) => x.norm() === 1n)).toBe(true);
    // closed under multiplication
    const names = new Set(u.map(s));
    expect(u.every((x) => u.every((y) => names.has(s(x.mul(y)))))).toBe(true);
  });
  test("division with a smaller remainder always works", () => {
    let seed = 12345;
    const rnd = () => {
      seed = (seed * 1103515245 + 12345) % 2147483648;
      return BigInt((seed >> 8) % 41) - 20n;
    };
    for (let t = 0; t < 400; t++) {
      const par = BigInt(t % 2);
      const mk = () => new C.Hurwitz([0, 0, 0, 0].map(() => 2n * rnd() + par));
      const x = mk();
      const y = mk();
      if (y.isZero()) continue;
      const { quotient, remainder } = x.divmod(y);
      expect(quotient.mul(y).add(remainder).equals(x)).toBe(true);
      expect(remainder.norm() < y.norm()).toBe(true);
    }
    expect(() => hurwitz(1, 0, 0, 0).divmod(hurwitz(0, 0, 0, 0))).toThrow(/zero/);
  });
  test("right greatest common divisors", () => {
    const c = hurwitz(1, 2, 0, 0);
    const g = C.hurwitzGcd(hurwitz(3, 1, -2, 5).mul(c), hurwitz(2, 1, 1, 0).mul(c));
    expect(g.norm()).toBe(5n);
    // g is c up to a unit on the left
    expect(C.hurwitzUnits().some((u) => u.mul(c).equals(g))).toBe(true);
  });
  test("a quaternion of prime norm", () => {
    for (const p of [2, 3, 5, 7, 13, 97, 1000003]) {
      const q = C.quaternionOfNorm(p);
      expect(q.isLipschitz).toBe(true);
      expect(q.norm()).toBe(BigInt(p));
    }
    expect(() => C.quaternionOfNorm(15)).toThrow(RangeError);
  });
  test("four squares", () => {
    for (let n = 0; n <= 400; n++) {
      const f = C.fourSquares(n);
      expect(f.reduce((t, x) => t + x * x, 0n)).toBe(BigInt(n));
    }
    expect(C.fourSquares(2026)).toEqual([23n, 23n, 22n, 22n]);
    const big = 10n ** 18n + 9n;
    expect(C.fourSquares(big).reduce((t, x) => t + x * x, 0n)).toBe(big);
    expect(() => C.fourSquares(-1)).toThrow(RangeError);
  });
  test("Jacobi's count", () => {
    const brute = (n) => {
      let c = 0;
      const m = Math.floor(Math.sqrt(n));
      for (let a = -m; a <= m; a++) for (let b = -m; b <= m; b++) for (let d = -m; d <= m; d++) {
        const r = n - a * a - b * b - d * d;
        if (r < 0) continue;
        const e = Math.round(Math.sqrt(r));
        if (e * e === r) c += e === 0 ? 1 : 2;
      }
      return BigInt(c);
    };
    for (const n of [1, 2, 7, 12, 30, 64]) expect(C.fourSquaresCount(n)).toBe(brute(n));
    expect(C.fourSquaresCount(0)).toBe(1n);
  });
});
