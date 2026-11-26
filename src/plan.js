// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Phases 1, 2 and 5 are built; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
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
    "phase": 6,
    "contents": "Complex gamma and zeta",
    "estimate": "700"
  }
].map(Object.freeze));
