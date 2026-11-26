# Changelog

## Unreleased

- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

## 0.1.0

Phases 1 and 2 of the roadmap.

- `complex`, `Complex`, `CQ`: exact complex rationals, read from text, with
  norms, exact moduli and exact square roots when they are rational.
- `rect`, `Rect`, `disc`, `Disc`: enclosures with rational data and arithmetic
  on them.
- `polar`, `modulus`, `argument`, `arg`, `exactTurn`, `fromPolar`, `roots`,
  `rootsOfUnity`.
- `exp sin cos tan sinh cosh tanh`: certified at any precision, on exact
  numbers, rectangles and discs.
- `log sqrt root pow asin acos atan asinh acosh atanh`: `Branched` results that
  report a branch cut, the value on each side and the conventional side;
  `{ branch, side }` options.
- `continueAlong`, `winding`: branches carried along polygonal paths, with
  exact crossing tests.
- `refined`: any function as a Core `Refinement`.
- `@ratmath/complex/cayley`: `cayleyDickson`, `quaternions`, `octonions`,
  `sedenions`, `Hyper`, with `properties()`; `rotate`, `rotationMatrix`,
  `fromRotationMatrix`, `fromAxis`, `rotationOf`, `primitive`, `interpolate`;
  `Hurwitz`, `hurwitz`, `hurwitzUnits`, `hurwitzGcd`, `quaternionOfNorm`,
  `fourSquares`, `fourSquaresCount`.
- Docs: eight topic pages, a reference and an explore page.
