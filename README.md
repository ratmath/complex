# @ratmath/complex

Complex numbers and complex analysis for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: set up, nothing built yet.** This repository holds the package
skeleton and the docs site. The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#3-a-complex-package-ratmathcomplex).

Exact complex rationals, elementary functions that say which branch they took, and complex analysis with proofs: contour integrals, residues and zero counts enclosed in rectangles that are guaranteed to contain the truth. Quaternions and octonions live here too.

## Planned

| Phase | Contents |
|---|---|
| 1 | Exact complex rationals, rectangles and discs, polar form, elementary functions with branch reporting |
| 2 | Cayley–Dickson numbers, quaternion rotations, integer quaternions |
| 3 | Complex Taylor models; contour integrals; argument principle and zero isolation |
| 4 | Residues, Laurent series, singularities, real integrals by residues |
| 5 | Möbius and conformal maps; graphics hooks |
| 6 | Complex gamma and zeta |

It will depend on: Core and analytic (real enclosures of exp, sin, atan and the constants), and optionally algebraic for conversions.

## Layout

- `index.js`, `src/`: the library. For now it exports only `plan`, the list
  above, so that the tests and the docs build have something to load.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, then one
  page per topic as topics are built). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/complex.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
