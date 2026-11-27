// @ratmath/complex: exact complex rationals, enclosures, and elementary
// functions that report their branches; Möbius maps, circles and the
// Riemann sphere; complex Taylor models, proven contour integrals, the
// argument principle and zeros with multiplicities; residues, Laurent series
// and real integrals by residues; gamma and zeta on the complex plane. Quaternions, octonions and the rest
// of the Cayley–Dickson numbers are in "@ratmath/complex/cayley".

export { Complex, complex, CQ, ComplexError } from "./src/complex.js";
export { Rect, rect, Disc, disc } from "./src/enclosure.js";
export {
  modulus, argument, arg, polar, fromPolar, exactTurn,
  exp, log, sqrt, root, roots, rootsOfUnity, pow,
  sin, cos, tan, sinh, cosh, tanh,
  asin, acos, atan, asinh, acosh, atanh,
  Branched, continueAlong, winding, refined,
} from "./src/functions.js";
export {
  Mobius, mobius, mobiusFrom, infinity, isInfinity, crossRatio, concyclic,
  GCircle, circle, circleThrough, lineThrough, toSphere, fromSphere, chordal2,
  cayley, discAutomorphism, halfPlaneToDisc, joukowski, regions, gridImage, sampleGrid,
} from "./src/mobius.js";
import { parse } from "@ratmath/analytic";
import { _setParser } from "./src/taylor.js";
_setParser(parse);
export { CTM, taylorModel } from "./src/taylor.js";
export { Contour, contourIntegral, residueSum, argumentPrinciple, countZeros, zeros } from "./src/contour.js";
export { laurent, residue, singularities, realLineIntegral, trigIntegral } from "./src/residues.js";
export { gamma, rgamma, logGamma, zeta, zetaTimesSMinus1, zetaZeros } from "./src/special.js";
export { cauchyRiemann } from "./src/cauchy-riemann.js";
export { plan } from "./src/plan.js";
