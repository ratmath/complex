# Changelog

## Unreleased

- `cauchyRiemann(u, v)`: checks the Cauchy–Riemann equations for f = u + iv
  given as formulas in x and y, with f′ when they hold and a certified
  witness point when they fail.

## 0.3.0

Phases 3, 4 and 6 of the roadmap, which completes it for this package.

- `taylorModel`, `CTM`: complex Taylor models, a polynomial with a proven bound
  on what it leaves out over a whole disc. Arithmetic, a checked reciprocal,
  `exp log sqrt root pow sin cos tan sinh cosh tanh` of a model, `at`, `range`,
  `derivative` (Cauchy's estimate), `integral`, and `zeroCount` (Rouché).
  The principal branches are refused on a disc that reaches their cut.
- `Contour`, `contourIntegral`, `residueSum`: paths of segments and circular
  arcs through exact points, integrated to a rectangle certain to contain the
  answer, with subdivision where the integrand needs it.
- `argumentPrinciple`, `countZeros`: zeros minus poles inside a contour as a
  proven whole number, from the winding of the image, with no derivative.
- `zeros`: every zero of an analytic function in a rectangle, isolated, with
  multiplicities, by exclusion, Rouché's theorem and an interval Newton step;
  the count is proved and what could not be settled is reported.
- `laurent`, `residue`, `singularities`: Laurent coefficients as rectangles,
  residues, and the kind of each singularity of a formula (pole with its
  order, removable, essential, branch point) from proven counts.
- `realLineIntegral`, `trigIntegral`: integrals of rational functions over the
  real line (with an optional e^(iax)) and of rational functions of cos and
  sin over a period, by residues, with the steps of the argument returned.
- `gamma`, `rgamma`, `logGamma`, `zeta`, `zetaTimesSMinus1`: the gamma and
  Riemann zeta functions for complex arguments at any precision, on points,
  enclosures and Taylor models. `zetaZeros` isolates the zeros in a range of
  heights and proves each to lie on the critical line.
- `Disc.toString` writes a centre with long fractions in decimals.
- `plan` is now empty.
- Docs: six topic pages (Taylor models, contour integrals, zeros, residues,
  real integrals, gamma and zeta) and an explore panel for zeros and poles.

- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

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
