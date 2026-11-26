# @ratmath/complex

Complex numbers and complex analysis for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

Exact complex rationals, certified enclosures of what is not rational, and
elementary functions that report their branch cuts instead of silently picking
a side. Quaternions, octonions and the rest of the Cayley–Dickson numbers are
at `@ratmath/complex/cayley`. Docs: [ratmath.com/complex](https://ratmath.com/complex/).

```js
import { complex, exp, log, sqrt } from "@ratmath/complex";

String(complex("(1 + 2i)/(3 - i)"));   // "1/10 + 7/10i", exactly
String(sqrt("3 + 4i"));                // "2 + i"
String(exp("1 + i"));                  // "1.4686939399 + 2.2873552872i ± 1.7e-21"
String(log("-1"));                     // on the cut: iπ from above, −iπ from below; convention: above
String(log("-1", { side: "below" }));  // "0 - 3.1415926536i ± 1.7e-21"

import { quaternion, rotate, fourSquares } from "@ratmath/complex/cayley";

String(quaternion("1 + 2i").mul("j"));    // "j + 2k"
rotate("1 + k", [1, 0, 0]).map(String);   // ["0", "1", "0"]
fourSquares(2026);                        // [23n, 23n, 22n, 22n]
```

## What is built

- **Complex rationals**: `complex`, `Complex`, the field object `CQ`. Exact
  arithmetic, norms, exact square roots when they exist.
- **Enclosures**: `rect`, `Rect`, `disc`, `Disc`, with rational data, and
  arithmetic that contains every possible result.
- **Polar form**: `polar`, `modulus`, `argument`, `exactTurn`, `fromPolar`,
  `roots`, `rootsOfUnity`.
- **Functions**: `exp sin cos tan sinh cosh tanh` at any precision, on exact
  numbers or enclosures.
- **Branches**: `log sqrt root pow asin acos atan asinh acosh atanh` return a
  `Branched` result that says when the input is on a cut and gives the value on
  each side; `continueAlong` and `winding` follow a branch along a path.
- **Cayley–Dickson numbers** (`/cayley`): `quaternions`, `octonions`,
  `sedenions`, `cayleyDickson(level, field)`, with `properties()` that computes
  which laws hold and exhibits zero divisors.
- **Rotations** (`/cayley`): `rotate`, `rotationMatrix`, `fromRotationMatrix`,
  `fromAxis`, `rotationOf`, `interpolate`, all exact.
- **Integer quaternions** (`/cayley`): `hurwitz`, `divmod`, `hurwitzGcd`,
  `quaternionOfNorm`, `fourSquares`, `fourSquaresCount`.

## Planned

From the [library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#3-a-complex-package-ratmathcomplex).

| Phase | Contents |
|---|---|
| 3 | Complex Taylor models; contour integrals; argument principle and zero isolation |
| 4 | Residues, Laurent series, singularities, real integrals by residues |
| 5 | Möbius and conformal maps; graphics hooks |
| 6 | Complex gamma and zeta |

It depends on Core and on analytic (certified real exp, log, sine, cosine and
arctangent). Analytic is not on npm yet, so it is a git submodule under
`vendor/` and a workspace: clone with `--recurse-submodules`.

## Layout

- `index.js`, `cayley.js`, `src/`: the library and its two entry points.
- `vendor/analytic`: the analytic library as a submodule.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, and one
  page per topic). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/complex.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
git submodule update --init
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
