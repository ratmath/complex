// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Nothing here is built yet; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 1,
    "contents": "Exact complex rationals, rectangles and discs, polar form, elementary functions with branch reporting",
    "estimate": "1,400"
  },
  {
    "phase": 2,
    "contents": "Cayley–Dickson numbers, quaternion rotations, integer quaternions",
    "estimate": "700"
  },
  {
    "phase": 3,
    "contents": "Complex Taylor models; contour integrals; argument principle and zero isolation",
    "estimate": "2,200"
  },
  {
    "phase": 4,
    "contents": "Residues, Laurent series, singularities, real integrals by residues",
    "estimate": "800"
  },
  {
    "phase": 5,
    "contents": "Möbius and conformal maps; graphics hooks",
    "estimate": "600, plus 800 in graphics"
  },
  {
    "phase": 6,
    "contents": "Complex gamma and zeta",
    "estimate": "700"
  }
].map(Object.freeze));
