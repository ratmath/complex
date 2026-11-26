// Phase 5: Möbius maps, generalized circles, the Riemann sphere, conformal maps.
// Exact answers are checked against floating-point eigenvalues, sampled
// points and identities computed another way.
import { describe, test, expect } from "bun:test";
import { Rational } from "@ratmath/core";
import {
  Complex, Mobius, mobius, mobiusFrom, infinity, crossRatio, concyclic, circle, circleThrough, lineThrough, GCircle,
  toSphere, fromSphere, chordal2, cayley, discAutomorphism, halfPlaneToDisc, joukowski, regions, gridImage, sampleGrid,
} from "../index.js";

const q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));
let seed = 7;
const rnd = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
const rq = () => q(Math.floor(rnd() * 21) - 10, Math.floor(rnd() * 5) + 1);
const rz = () => new Complex(rq(), rq());
const randomMap = () => {
  for (;;) {
    try {
      return new Mobius(rz(), rz(), rz(), rz());
    } catch {}
  }
};
// rational points on the circle |z − c|² = r² with r rational: c + r·(1 − t² + 2ti)/(1 + t²)
const onCircle = (c, r, t) => Complex.from(c).add(new Complex(q(1).subtract(t.multiply(t)), q(2).multiply(t)).scale(r.divide(q(1).add(t.multiply(t)))));
const fl = (z) => (z === infinity ? null : z.toNumbers());

describe("Möbius maps", () => {
  test("three points to three, composition and inverse", () => {
    for (let i = 0; i < 20; i++) {
      const from = [rz(), rz(), rz()];
      const to = [rz(), rz(), rz()];
      if (new Set(from.map(String)).size < 3 || new Set(to.map(String)).size < 3) continue;
      const T = mobiusFrom(from, to);
      from.forEach((z, k) => expect(T.apply(z).equals(to[k])).toBe(true));
      expect(T.compose(T.inverse()).isIdentity).toBe(true);
      // composition agrees with applying one after the other (floats as an independent route)
      const S = randomMap();
      const z = rz();
      const a = S.compose(T).apply(z);
      const b = S.applyFloat(T.applyFloat(z.toNumbers()) ?? [NaN, NaN]);
      if (a !== infinity && b) expect(Math.hypot(a.toNumbers()[0] - b[0], a.toNumbers()[1] - b[1])).toBeLessThan(1e-9 * (1 + Math.hypot(...b)));
    }
    const T = mobiusFrom([0, 1, infinity], ["i", 1, -1]);
    expect(T.apply(infinity).equals(Complex.from(-1))).toBe(true);
    expect(mobius(0, 1, 1, 0).apply(0)).toBe(infinity);
    expect(() => mobius(1, 2, 2, 4)).toThrow();
  });
  test("the cross-ratio is kept, and is real exactly for concyclic points", () => {
    for (let i = 0; i < 20; i++) {
      const zs = [rz(), rz(), rz(), rz()];
      if (new Set(zs.map(String)).size < 4) continue;
      const T = randomMap();
      const w = zs.map((z) => T.apply(z));
      const a = crossRatio(...zs);
      const b = crossRatio(...w);
      expect(a === infinity ? b === infinity : a.equals(b)).toBe(true);
    }
    const c = new Complex(q(1, 3), q(-2));
    const r = q(5, 7);
    const pts = [q(0), q(1), q(2, 3), q(-5)].map((t) => onCircle(c, r, t));
    expect(concyclic(...pts)).toBe(true);
    expect(concyclic(pts[0], pts[1], pts[2], c)).toBe(false);
    expect(concyclic(0, 1, 2, infinity)).toBe(true); // a line through ∞
  });
  test("the kind of map agrees with the eigenvalues", () => {
    const kindFromEigen = (T) => {
      // eigenvalues of [[a, b], [c, d]] in floats, k = λ₁/λ₂
      const [a, b, c, d] = [T.a, T.b, T.c, T.d].map((z) => z.toNumbers());
      const tr = [a[0] + d[0], a[1] + d[1]];
      const det = T.det.toNumbers();
      const disc = [tr[0] * tr[0] - tr[1] * tr[1] - 4 * det[0], 2 * tr[0] * tr[1] - 4 * det[1]];
      const m = Math.hypot(...disc);
      const s = [Math.sqrt((m + disc[0]) / 2), Math.sign(disc[1] || 1) * Math.sqrt(Math.max(0, (m - disc[0]) / 2))];
      const l1 = [(tr[0] + s[0]) / 2, (tr[1] + s[1]) / 2];
      const l2 = [(tr[0] - s[0]) / 2, (tr[1] - s[1]) / 2];
      const n2 = l2[0] ** 2 + l2[1] ** 2;
      const k = [(l1[0] * l2[0] + l1[1] * l2[1]) / n2, (l1[1] * l2[0] - l1[0] * l2[1]) / n2];
      if (Math.hypot(k[0] - 1, k[1]) < 1e-9) return "parabolic";
      if (Math.abs(Math.hypot(...k) - 1) < 1e-9) return "elliptic";
      if (Math.abs(k[1]) < 1e-9 && k[0] > 0) return "hyperbolic";
      return "loxodromic";
    };

    for (let i = 0; i < 60; i++) {
      const T = randomMap();
      if (T.isIdentity) continue;
      expect(T.classify()).toBe(kindFromEigen(T));

    }
    // conjugates of z ↦ kz and z ↦ z + 1 for each kind
    const P = randomMap();
    for (const [k, kind] of [["3/5 + 4/5i", "elliptic"], ["4", "hyperbolic"], ["-2", "loxodromic"], ["1 + i", "loxodromic"]]) {
      const T = P.compose(mobius(k, 0, 0, 1)).compose(P.inverse());
      expect(T.classify()).toBe(kind);
    }
    expect(P.compose(mobius(1, 1, 0, 1)).compose(P.inverse()).classify()).toBe("parabolic");
  });
  test("fixed points", () => {
    for (let i = 0; i < 30; i++) {
      const T = randomMap();
      if (T.isIdentity) continue;
      const f = T.fixedPoints();
      if (f.exact) {
        for (const p of f.points) {
          const w = T.apply(p);
          expect(p === infinity ? w === infinity : w !== infinity && w.equals(p)).toBe(true);
        }
      } else {
        for (const R of f.points) {
          const c = R.center.toNumbers();
          const w = T.applyFloat(c);
          expect(Math.hypot(w[0] - c[0], w[1] - c[1])).toBeLessThan(1e-9 * (1 + Math.hypot(...c)));
        }
      }
    }
    expect(mobius(2, 0, 0, 1).fixedPoints().points.map(String)).toEqual(["0", "∞"]);
    expect(mobius(1, 1, 0, 1).fixedPoints().points).toEqual([infinity]);
    const g = mobius(1, 1, 1, 0).fixedPoints(); // z = 1 + 1/z: the golden ratio
    expect(g.exact).toBe(false);
    expect(g.points[0].contains(Complex.from(0)) ).toBe(false);
    expect(Math.abs(g.points[0].center.toNumbers()[0] - (1 + Math.sqrt(5)) / 2)).toBeLessThan(1e-15);
  });
});

describe("generalized circles", () => {
  test("images of circles contain the images of their points", () => {
    for (let i = 0; i < 15; i++) {
      const c = rz();
      const r = q(Math.floor(rnd() * 9) + 1, Math.floor(rnd() * 4) + 1);
      const K = circle(c, r.multiply(r));
      const T = randomMap();
      const img = K.map(T);
      for (const t of [q(0), q(1), q(-3, 2), q(7, 5)]) {
        const p = onCircle(c, r, t);
        expect(K.contains(p)).toBe(true);
        expect(img.contains(T.apply(p))).toBe(true);
      }
    }
  });
  test("lines, three points, and the Cayley map", () => {
    expect(circleThrough(1, "i", -1).equals(circle(0, 1))).toBe(true);
    expect(circleThrough(0, 1, 2).isLine).toBe(true);
    // Cayley sends the real line to the unit circle and the upper half-plane into the disc
    const C = cayley();
    expect(lineThrough(0, 1).map(C).equals(circle(0, 1))).toBe(true);
    for (let i = 0; i < 20; i++) {
      const z = new Complex(rq(), q(Math.floor(rnd() * 9) + 1, 3));
      expect(regions.disc(C.apply(z))).toBe(true);
    }
    expect(C.apply("i").isZero()).toBe(true);
  });
  test("grid images are exact", () => {
    const T = randomMap();
    const g = gridImage(T, { re: [-1, 1], im: [-1, 1], step: "1/2" });
    expect(g.lines.length).toBe(10);
    for (const { from, image } of g.lines) {
      const v = q(...from.split(" = ")[1].split("/").map(Number));
      const pts = [q(-2), q(1, 3), q(5)].map((t) => (from.startsWith("re") ? new Complex(v, t) : new Complex(t, v)));
      for (const p of pts) {
        const w = T.apply(p);
        expect(w === infinity ? image.isLine : image.contains(w)).toBe(true);
      }
    }
  });
});

describe("the Riemann sphere and conformal maps", () => {
  test("stereographic projection is exact and invertible", () => {
    for (let i = 0; i < 20; i++) {
      const z = rz();
      const P = toSphere(z);
      expect(P.reduce((s, x) => s.add(x.multiply(x)), q(0)).equals(q(1))).toBe(true);
      expect(fromSphere(P).equals(z)).toBe(true);
    }
    expect(fromSphere([0, 0, 1])).toBe(infinity);
    expect(chordal2(0, infinity).equals(q(4))).toBe(true);
  });
  test("rotations of the sphere keep chordal distances", () => {
    // unitary up to scale: [[a, b], [−b̄, ā]]
    for (let i = 0; i < 10; i++) {
      const a = rz();
      const b = rz();
      if (a.isZero() && b.isZero()) continue;
      const R = mobius(a, b, b.conj().neg(), a.conj());
      expect(R.isSphereRotation).toBe(true);
      const z = rz();
      const w = rz();
      const [Rz, Rw] = [R.apply(z), R.apply(w)];
      expect(chordal2(Rz, Rw).equals(chordal2(z, w))).toBe(true);
    }
    expect(mobius(2, 0, 0, 1).isSphereRotation).toBe(false);
  });
  test("disc automorphisms keep the disc and its boundary", () => {
    for (let i = 0; i < 10; i++) {
      const a = new Complex(q(Math.floor(rnd() * 9) - 4, 10), q(Math.floor(rnd() * 9) - 4, 10));
      const D = discAutomorphism(a, "3/5 + 4/5i");
      expect(D.apply(a).isZero()).toBe(true);
      for (const t of [q(0), q(2), q(-1, 3)]) expect(D.apply(onCircle(0, q(1), t)).norm().equals(q(1))).toBe(true);
      expect(regions.disc(D.apply(new Complex(q(1, 3), q(-1, 4))))).toBe(true);
    }
    expect(() => discAutomorphism(2)).toThrow();
    expect(() => discAutomorphism(0, 2)).toThrow();
    const H = halfPlaneToDisc("1 + 2i");
    expect(H.apply("1 + 2i").isZero()).toBe(true);
    expect(H.apply(5).norm().equals(q(1))).toBe(true);
  });
  test("Joukowski sends circles |z| = r to ellipses with semi-axes (r ± 1/r)/2", () => {
    const r = q(2);
    for (const t of [q(0), q(1), q(1, 2), q(-3)]) {
      const w = joukowski(onCircle(0, r, t));
      const A = r.add(q(1).divide(r)).divide(q(2));
      const B = r.subtract(q(1).divide(r)).divide(q(2));
      expect(w.re.multiply(w.re).divide(A.multiply(A)).add(w.im.multiply(w.im).divide(B.multiply(B))).equals(q(1))).toBe(true);
    }
    expect(joukowski(0)).toBe(infinity);
  });
  test("sampled grids agree with the exact map", () => {
    const T = randomMap();
    const s = sampleGrid(T, { width: 8, height: 8 });
    expect(s.kind).toBe("estimated");
    const z = new Complex(q(-2).add(q(1, 4)), q(-2).add(q(1, 4))); // the centre of the first pixel
    const w = T.apply(z);
    if (w !== infinity) expect(Math.hypot(s.values[0] - w.toNumbers()[0], s.values[1] - w.toNumbers()[1])).toBeLessThan(1e-9 * (1 + Math.hypot(...w.toNumbers())));
    expect(sampleGrid(([x, y]) => [x * x - y * y, 2 * x * y], { width: 2, height: 2 }).values.length).toBe(8);
  });
});

test("fixed points when the discriminant is a negative real (on the square root's cut)", () => {
  // disc = −12: fixed points ±√3 i·… must still be tight
  const T = mobiusFrom([0, 1, infinity], [-1, "i", 1]);
  const f = T.fixedPoints();
  expect(f.exact).toBe(false);
  for (const R of f.points) {
    expect(R.width.toNumber()).toBeLessThan(1e-15);
    const c = R.center.toNumbers();
    const w = T.applyFloat(c);
    expect(Math.hypot(w[0] - c[0], w[1] - c[1])).toBeLessThan(1e-12);
  }
});
