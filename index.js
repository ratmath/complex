// @ratmath/complex: exact complex rationals, enclosures, and elementary
// functions that report their branches. Quaternions, octonions and the rest
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
export { plan } from "./src/plan.js";
