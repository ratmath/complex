# Changelog

## Unreleased

## 0.2.0

Phase 5 of the roadmap.

- `mobius`, `Mobius`, `mobiusFrom`, `crossRatio`, `concyclic`, `infinity`:
  exact Möbius maps of the Riemann sphere, classified by (a + d)²/(ad − bc),
  with fixed points exact when the square root is rational and certified
  otherwise.
- `GCircle`, `circle`, `circleThrough`, `lineThrough`: generalized circles as
  Hermitian forms with rational data, mapped exactly.
- `toSphere`, `fromSphere`, `chordal2`: stereographic projection, exact.
- `cayley`, `discAutomorphism`, `halfPlaneToDisc`, `joukowski`, `regions`.
- `gridImage` (exact circles) and `sampleGrid` (estimated values): data for
  pictures in graphics.
- Docs: Möbius maps and Riemann sphere topic pages, an explore panel.

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
