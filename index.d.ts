// Type declarations for @ratmath/complex.

import type { Rational, RationalInterval } from "@ratmath/core";
import type { Refinement } from "@ratmath/core/refine";

/** Anything read as a real rational: a Rational, bigint, number, or text in Core's notation. */
export type RationalLike = Rational | bigint | number | string;
/** Anything read as an exact complex number. */
export type ComplexLike = Complex | RationalLike | readonly [RationalLike, RationalLike] | { re: RationalLike; im: RationalLike };
/** An exact number or an enclosure. */
export type Enclosable = ComplexLike | Rect | Disc;

export type ComplexErrorCode = "division-by-zero" | "branch-cut" | "branch-point" | "pole" | "precision" | "bad-option" | "no-quotient" | "not-a-rotation" | "too-large" | "bad-argument" | "degenerate" | "singularity" | "zero-on-contour" | "unsupported";
export class ComplexError extends Error {
  code: ComplexErrorCode;
  constructor(code: ComplexErrorCode, message: string);
}

export class Complex {
  readonly re: Rational;
  readonly im: Rational;
  constructor(re?: RationalLike, im?: RationalLike);
  static from(x: ComplexLike): Complex;
  /** Exact arithmetic on numbers, i, + - * /, ^ with a whole-number exponent, and parentheses. */
  static parse(text: string): Complex;
  static readonly ZERO: Complex;
  static readonly ONE: Complex;
  static readonly I: Complex;
  add(o: ComplexLike): Complex;
  sub(o: ComplexLike): Complex;
  mul(o: ComplexLike): Complex;
  div(o: ComplexLike): Complex;
  neg(): Complex;
  conj(): Complex;
  inv(): Complex;
  /** |z|², exact. */
  norm(): Rational;
  /** |z| when it is rational, else null. */
  exactModulus(): Rational | null;
  pow(n: number | bigint): Complex;
  scale(q: RationalLike): Complex;
  equals(o: ComplexLike): boolean;
  isZero(): boolean;
  isReal(): boolean;
  isImaginary(): boolean;
  isGaussianInteger(): boolean;
  /** [w, −w] with the principal root first when the roots are complex rationals, else null. */
  exactSqrt(): [Complex, Complex] | null;
  toNumbers(): [number, number];
  toString(): string;
  toTeX(): string;
  toJSON(): { re: string; im: string };
}
export function complex(re: ComplexLike, im?: RationalLike): Complex;

/** A RatMath field object. */
export interface Field<T> {
  name?: string;
  zero: T;
  one: T;
  add(a: T, b: T): T;
  sub(a: T, b: T): T;
  mul(a: T, b: T): T;
  div(a: T, b: T): T;
  neg(a: T): T;
  sign(a: T): number;
  from?(x: unknown): T;
  format?(a: T): string;
}
/** The complex rationals Q(i) as a field object. */
export const CQ: Field<Complex> & { inv(a: ComplexLike): Complex; equals(a: ComplexLike, b: ComplexLike): boolean; of(x: ComplexLike): Complex; conj(a: ComplexLike): Complex; characteristic: 0 };

export interface Bounds {
  lo: Rational;
  hi: Rational;
}
export type IntervalLike = RationalLike | RationalInterval | Bounds | readonly [RationalLike, RationalLike];

export class Rect {
  readonly re: Bounds;
  readonly im: Bounds;
  constructor(re: Bounds, im: Bounds);
  static from(x: Enclosable): Rect;
  static hull(...items: Enclosable[]): Rect;
  readonly isPoint: boolean;
  readonly center: Complex;
  readonly width: Rational;
  readonly height: Rational;
  readonly size: Rational;
  readonly real: RationalInterval;
  readonly imag: RationalInterval;
  corners(): Complex[];
  add(o: Enclosable): Rect;
  sub(o: Enclosable): Rect;
  mul(o: Enclosable): Rect;
  div(o: Enclosable): Rect;
  inv(): Rect;
  neg(): Rect;
  conj(): Rect;
  scale(q: RationalLike): Rect;
  timesI(): Rect;
  pow(n: number): Rect;
  normRange(): Bounds;
  containsZero(): boolean;
  contains(o: Enclosable): boolean;
  intersects(o: Enclosable): boolean;
  intersect(o: Enclosable): Rect | null;
  hull(o: Enclosable): Rect;
  round(bits?: number): Rect;
  toDisc(): Disc;
  split(): Rect[];
  toString(digits?: number): string;
  toJSON(): { re: [string, string]; im: [string, string] };
}
/** rect("1 + 2i"), rect([0, 1], [2, 3]), or rect(z, radius) for the square around z. */
export function rect(a: Enclosable | IntervalLike, b?: IntervalLike): Rect;

export class Disc {
  readonly center: Complex;
  readonly radius: Rational;
  constructor(center: ComplexLike, radius?: RationalLike);
  static from(x: Enclosable): Disc;
  readonly isPoint: boolean;
  add(o: Enclosable): Disc;
  sub(o: Enclosable): Disc;
  mul(o: Enclosable): Disc;
  div(o: Enclosable): Disc;
  inv(): Disc;
  neg(): Disc;
  conj(): Disc;
  scale(q: RationalLike): Disc;
  contains(o: Enclosable): boolean;
  containsZero(): boolean;
  intersects(o: Enclosable): boolean;
  hull(o: Enclosable): Disc;
  toRect(): Rect;
  toString(digits?: number): string;
}
export function disc(center: ComplexLike, radius?: RationalLike): Disc;

export type Side = "above" | "below" | "right" | "left";
export interface BranchOptions {
  bits?: number;
  /** A numbered branch: adds 2πik to a logarithm, multiplies an n-th root by e^(2πik/n). */
  branch?: number;
  /** Which side of the cut the input is on. */
  side?: Side;
}
export type Precision = number | { bits?: number };
export interface Candidate {
  side: Side;
  value: Rect;
}
/** The result of a function with a branch cut. */
export class Branched {
  readonly name: string;
  readonly cut: string;
  readonly value: Rect | null;
  readonly straddles: boolean;
  readonly candidates: readonly Candidate[];
  readonly convention: Candidate | null;
  readonly branch: "principal" | number;
  /** The enclosure, or a ComplexError with code "branch-cut". */
  readonly rect: Rect;
  /** A rectangle around every candidate. */
  readonly hull: Rect;
  toString(digits?: number): string;
}

export function modulus(z: Enclosable, bits?: number): RationalInterval;
export function argument(z: ComplexLike, bits?: number): RationalInterval;
export function arg(z: Enclosable, opts?: number | BranchOptions): Branched;
export function exactTurn(z: ComplexLike): Rational | null;
export function polar(z: ComplexLike, bits?: number): { modulus: RationalInterval; argument: RationalInterval; normSquared: Rational; exactModulus: Rational | null; exactTurn: Rational | null };
export function fromPolar(r: IntervalLike, theta: IntervalLike | { turns: RationalLike }, bits?: number): Rect;

export function exp(z: Enclosable, opts?: Precision): Rect;
export function exp(z: CTM): CTM;
export function sin(z: Enclosable, opts?: Precision): Rect;
export function sin(z: CTM): CTM;
export function cos(z: Enclosable, opts?: Precision): Rect;
export function cos(z: CTM): CTM;
export function tan(z: Enclosable, opts?: Precision): Rect;
export function tan(z: CTM): CTM;
export function sinh(z: Enclosable, opts?: Precision): Rect;
export function sinh(z: CTM): CTM;
export function cosh(z: Enclosable, opts?: Precision): Rect;
export function cosh(z: CTM): CTM;
export function tanh(z: Enclosable, opts?: Precision): Rect;
export function tanh(z: CTM): CTM;

export function log(z: Enclosable, opts?: number | BranchOptions): Branched;
/** Of a Taylor model: the model of the principal branch (or { branch: k }); refused when the disc reaches the cut. */
export function log(z: CTM, opts?: { branch?: number }): CTM;
export function sqrt(z: Enclosable, opts?: number | BranchOptions): Branched;
export function sqrt(z: CTM, opts?: { branch?: number }): CTM;
export function root(z: Enclosable, n: number, opts?: number | BranchOptions): Branched;
export function root(z: CTM, n: number, opts?: { branch?: number }): CTM;
export function pow(z: Enclosable, w: Enclosable, opts?: number | BranchOptions): Branched;
export function pow(z: CTM, w: Enclosable | CTM): CTM;
export function pow(z: Enclosable, w: CTM): CTM;
export function roots(z: ComplexLike, n: number, opts?: Precision): Rect[];
export function rootsOfUnity(n: number, bits?: number): Rect[];
export function asin(z: Enclosable, opts?: number | BranchOptions): Branched;
export function acos(z: Enclosable, opts?: number | BranchOptions): Branched;
export function atan(z: Enclosable, opts?: number | BranchOptions): Branched;
export function asinh(z: Enclosable, opts?: number | BranchOptions): Branched;
export function acosh(z: Enclosable, opts?: number | BranchOptions): Branched;
export function atanh(z: Enclosable, opts?: number | BranchOptions): Branched;

/** The winding number of a closed polygon about w (default 0), exact. */
export function winding(path: readonly ComplexLike[], w?: ComplexLike): number;
/** Carry a branch of log, a root or a power continuously along a polygonal path. */
export function continueAlong(
  fn: "log" | { root: number } | { pow: ComplexLike },
  path: readonly ComplexLike[],
  opts?: { bits?: number; branch?: number; about?: ComplexLike },
): { value: Rect; turns: number; branch: number; values: Rect[] };

export type FunctionName = "exp" | "sin" | "cos" | "tan" | "sinh" | "cosh" | "tanh" | "log" | "sqrt" | "asin" | "acos" | "atan" | "asinh" | "acosh" | "atanh";
/** A function value as a Core Refinement: a quick preview, then the certified enclosure. */
export function refined(name: FunctionName, z: Enclosable, opts?: BranchOptions): Refinement;

// ---- Möbius maps, generalized circles, the Riemann sphere ----

/** The point at infinity of the Riemann sphere. */
export const infinity: { readonly infinity: true; toString(): string; toJSON(): string };
export type Infinity = typeof infinity;
/** A point of the Riemann sphere: a complex number or ∞. */
export type SpherePoint = ComplexLike | Infinity;
export function isInfinity(z: unknown): boolean;

export type MobiusKind = "identity" | "elliptic" | "parabolic" | "hyperbolic" | "loxodromic";
/** z ↦ (az + b)/(cz + d) with complex rational coefficients and ad − bc ≠ 0. Exact. */
export class Mobius {
  readonly a: Complex;
  readonly b: Complex;
  readonly c: Complex;
  readonly d: Complex;
  constructor(a: ComplexLike, b: ComplexLike, c: ComplexLike, d: ComplexLike);
  readonly det: Complex;
  readonly matrix: readonly [readonly [Complex, Complex], readonly [Complex, Complex]];
  apply(z: SpherePoint): Complex | Infinity;
  compose(o: Mobius): Mobius;
  inverse(): Mobius;
  /** Equal as maps (coefficients up to a common factor). */
  equals(o: Mobius): boolean;
  /** (a + d)²/(ad − bc), the invariant that classifies the map. */
  readonly traceSquared: Complex;
  readonly isIdentity: boolean;
  classify(): MobiusKind;
  fixedPoints(opts?: { bits?: number }): { exact: true; points: (Complex | Infinity)[] | "every point" } | { exact: false; points: Rect[] };
  /** Is it a rotation of the Riemann sphere (unitary up to scale)? */
  readonly isSphereRotation: boolean;
  /** Estimated: the map on a pair of doubles, null at the pole. */
  applyFloat(p: readonly [number, number]): [number, number] | null;
  toString(): string;
}
export function mobius(a: ComplexLike, b: ComplexLike, c: ComplexLike, d: ComplexLike): Mobius;
/** The unique Möbius map sending three distinct points to three distinct points. */
export function mobiusFrom(from: readonly [SpherePoint, SpherePoint, SpherePoint], to: readonly [SpherePoint, SpherePoint, SpherePoint]): Mobius;
export function crossRatio(z1: SpherePoint, z2: SpherePoint, z3: SpherePoint, z4: SpherePoint): Complex | Infinity;
/** Do four points lie on one circle or line? Exact. */
export function concyclic(z1: SpherePoint, z2: SpherePoint, z3: SpherePoint, z4: SpherePoint): boolean;

/** A circle or line, as the Hermitian form A|z|² + B̄z + Bz̄ + C = 0 (A, C real). */
export class GCircle {
  readonly A: Rational;
  readonly B: Complex;
  readonly C: Rational;
  constructor(A: RationalLike, B: ComplexLike, C: RationalLike);
  readonly isLine: boolean;
  readonly center: Complex | null;
  readonly radius2: Rational | null;
  value(z: ComplexLike): Rational;
  contains(z: SpherePoint): boolean;
  equals(o: GCircle): boolean;
  /** The image under a Möbius map, exactly. */
  map(T: Mobius): GCircle;
  toString(): string;
}
export function circle(center: ComplexLike, radius2: RationalLike): GCircle;
export function circleThrough(z1: ComplexLike, z2: ComplexLike, z3: ComplexLike): GCircle;
export function lineThrough(z1: ComplexLike, z2: ComplexLike): GCircle;

/** Stereographic projection to the unit sphere, exact. */
export function toSphere(z: SpherePoint): [Rational, Rational, Rational];
export function fromSphere(p: readonly [RationalLike, RationalLike, RationalLike]): Complex | Infinity;
/** The squared chordal distance on the Riemann sphere. */
export function chordal2(z: SpherePoint, w: SpherePoint): Rational;

/** The Cayley map z ↦ (z − i)/(z + i), upper half-plane onto the disc. */
export function cayley(): Mobius;
export function discAutomorphism(a: ComplexLike, u?: ComplexLike): Mobius;
export function halfPlaneToDisc(z0: ComplexLike): Mobius;
export function joukowski(z: ComplexLike): Complex | Infinity;
export const regions: {
  readonly disc: (z: SpherePoint) => boolean;
  readonly upperHalfPlane: (z: SpherePoint) => boolean;
  readonly rightHalfPlane: (z: SpherePoint) => boolean;
  readonly firstQuadrant: (z: SpherePoint) => boolean;
};
export interface GridRange { re?: readonly [RationalLike, RationalLike]; im?: readonly [RationalLike, RationalLike] }
/** Exact: the image of a grid of lines, each a generalized circle. */
export function gridImage(T: Mobius, opts?: GridRange & { step?: RationalLike }): { readonly kind: "exact"; readonly lines: readonly { from: string; image: GCircle }[] };
/** Estimated: values of a map on a grid of doubles. */
export function sampleGrid(
  f: Mobius | ((p: [number, number]) => [number, number] | null),
  opts?: { re?: readonly [number, number]; im?: readonly [number, number]; width?: number; height?: number },
): { readonly kind: "estimated"; readonly width: number; readonly height: number; readonly re: [number, number]; readonly im: [number, number]; readonly values: Float64Array };

// ---- Taylor models ----

/** A function of z: text ("exp(z) - 3z^2"), a function of one Taylor model, or a syntax tree from analytic's parser. */
export type AnalyticFunction = string | ((z: CTM) => CTM | Enclosable) | { type: string };
export interface ModelOptions {
  /** The degree kept (default 12). */
  order?: number;
  /** The working precision of the coefficients in binary digits (default 96). */
  bits?: number;
}
/**
 * A complex Taylor model: a polynomial P and a bound ε with |f(z) − P(z − c)| ≤ ε
 * on the whole closed disc |z − c| ≤ r, where f is analytic.
 */
export class CTM {
  readonly center: Complex;
  readonly radius: Rational;
  readonly order: number;
  readonly bits: number;
  readonly disc: Disc;
  /** The proven bound ε. */
  readonly error: Rational;
  /** The polynomial's coefficients in powers of (z − c). */
  readonly coefficients: Complex[];
  /** The model of z itself on a disc. */
  static z(disc: Disc | Rect | ComplexLike, opts?: ModelOptions): CTM;
  static constant(disc: Disc | Rect | ComplexLike, value: Enclosable, opts?: ModelOptions): CTM;
  add(o: CTM | Enclosable): CTM;
  sub(o: CTM | Enclosable): CTM;
  mul(o: CTM | Enclosable): CTM;
  div(o: CTM | Enclosable): CTM;
  neg(): CTM;
  /** 1/f; refuses ("division-by-zero") when f may vanish on the disc. */
  inv(): CTM;
  pow(n: number | bigint): CTM;
  scale(q: RationalLike): CTM;
  /** A rectangle certain to contain f(z), for a point of the disc. */
  at(z: ComplexLike): Rect;
  /** A disc certain to contain every value of f on the model's disc. */
  range(): Disc;
  /** The number of zeros in the disc when Rouché's theorem decides it from the coefficients, else null. */
  zeroCount(): number | null;
  /** A model of f′ on a strictly smaller disc (default half the radius). */
  derivative(radius?: RationalLike): CTM;
  /** A rectangle certain to contain the integral of f from a to b, two points of the disc. */
  integral(a: ComplexLike, b: ComplexLike): Rect;
  toString(digits?: number): string;
}
export function taylorModel(f: AnalyticFunction, disc: Disc | Rect | ComplexLike, opts?: ModelOptions): CTM;

// ---- contours, integrals, zeros ----

export interface ContourPiece { readonly a: Complex; readonly b: Complex; readonly center: Complex | null; readonly ccw: boolean }
export type ContourLike = Contour | Rect | Disc | readonly ComplexLike[];
export class Contour {
  readonly pieces: readonly ContourPiece[];
  readonly start: Complex;
  readonly end: Complex;
  readonly closed: boolean;
  static from(x: ContourLike): Contour;
  static path(points: readonly ComplexLike[]): Contour;
  static polygon(points: readonly ComplexLike[]): Contour;
  static rectangle(r: Rect | ComplexLike): Contour;
  static circle(center: ComplexLike, radius: RationalLike, opts?: { clockwise?: boolean }): Contour;
  static arc(center: ComplexLike, from: ComplexLike, to: ComplexLike, opts?: { clockwise?: boolean }): Contour;
  join(other: ContourLike): Contour;
  reversed(): Contour;
  /** Points along the contour as doubles, for drawing. */
  sample(perPiece?: number): [number, number][];
  /** How many times the closed contour goes round w. */
  winding(w?: ComplexLike): number;
  toString(): string;
}
export interface IntegralOptions extends ModelOptions {
  /** The width to aim for (default 2^−40). */
  tolerance?: RationalLike;
  maxDepth?: number;
  maxPieces?: number;
}
/** The integral of f along a contour, as a rectangle certain to contain it. */
export function contourIntegral(f: AnalyticFunction, contour: ContourLike, opts?: IntegralOptions): Rect;
/** The integral round a closed contour over 2πi: the sum of the residues inside. */
export function residueSum(f: AnalyticFunction, contour: ContourLike, opts?: IntegralOptions): Rect;
/** The winding number of f(contour) round 0: zeros minus poles inside, a proven whole number. */
export function argumentPrinciple(f: AnalyticFunction, contour: ContourLike, opts?: IntegralOptions): number;
export function countZeros(f: AnalyticFunction, contour: ContourLike, opts?: IntegralOptions): number;
export interface FoundZero {
  enclosure: Rect;
  /** The total multiplicity of the zeros in the enclosure. */
  multiplicity: number;
  /** True for one simple zero. */
  simple: boolean;
  /** A disc in which the zeros of the enclosure are the only ones. */
  isolating: Disc | null;
  /** The enclosure straddles the region's edge: listed, not counted. */
  boundary: boolean;
  possibleDuplicate?: boolean;
}
export interface Unresolved { region: Rect; reason: string }
export interface ZerosResult { zeros: FoundZero[]; count: number; complete: boolean; unresolved: Unresolved[] }
export interface ZerosOptions {
  /** Binary digits of each simple zero (default 64). */
  bits?: number;
  /** How closely two zeros are told apart (default 2^−20). */
  tolerance?: RationalLike;
  order?: number;
  modelBits?: number;
  maxCells?: number;
}
/** The zeros of an analytic f in a rectangle, isolated and counted with multiplicity. */
export function zeros(f: AnalyticFunction, region: Rect | ComplexLike, opts?: ZerosOptions): ZerosResult;

// ---- residues ----

export interface LaurentSeries {
  readonly center: Complex;
  readonly radius: Rational;
  readonly from: number;
  readonly to: number;
  readonly coefficients: readonly Rect[];
  coefficient(n: number): Rect;
  readonly principalPart: { n: number; value: Rect; nonzero: boolean }[];
  toString(digits?: number): string;
}
export function laurent(f: AnalyticFunction, center: ComplexLike, radius: RationalLike, opts?: IntegralOptions & { from?: number; to?: number }): LaurentSeries;
export function residue(f: AnalyticFunction, point: Enclosable, opts?: IntegralOptions & { radius?: RationalLike }): Rect;
export interface Singularity {
  enclosure: Rect;
  kind: "pole" | "removable" | "essential" | "branch point" | "undetermined";
  /** A pole's order; for a removable singularity the order of the zero left over; otherwise null. */
  order: number | null;
  note: string;
  boundary?: boolean;
}
/** The singularities of a formula (text) in a rectangle, with their kinds. */
export function singularities(f: string | { type: string }, region: Rect | ComplexLike, opts?: ZerosOptions): { singularities: Singularity[]; complete: boolean; unresolved: Unresolved[]; numerator: string; denominator: string };
export interface Step { title: string; detail: string }
export interface PoleRecord { enclosure: Rect; order: number; taken: boolean; residue: Rect | null }
export interface Bounds { lo: Rational; hi: Rational }
export interface RealIntegral {
  readonly value: Rect;
  /** The interval for a real integrand. */
  readonly real?: Bounds;
  /** With `fourier`: ∫ f cos(ax) and ∫ f sin(ax), for f real on the real axis. */
  readonly cos?: Bounds;
  readonly sin?: Bounds;
  readonly poles: PoleRecord[];
  readonly steps: Step[];
}
/** ∫ P(x)/Q(x) dx over the real line (times e^(iax) with `fourier: a`) by residues, with the steps. */
export function realLineIntegral(f: string | { type: string }, opts?: { fourier?: RationalLike; bits?: number; variable?: string; tolerance?: RationalLike }): RealIntegral;
/** ∫₀^2π of a rational function of cos t and sin t by residues, with the steps. */
export function trigIntegral(f: string | { type: string }, opts?: { bits?: number; variable?: string; tolerance?: RationalLike }): RealIntegral;

// ---- gamma and zeta ----

/** Γ(z). A point or an enclosure gives a Rect; a Taylor model gives a model. Refuses ("pole") at 0, −1, −2, … */
export function gamma(z: Enclosable, opts?: number | { bits?: number }): Rect;
export function gamma(z: CTM): CTM;
/** 1/Γ(z), entire. */
export function rgamma(z: Enclosable, opts?: number | { bits?: number }): Rect;
export function rgamma(z: CTM): CTM;
/** log Γ(z) for Re z > 0, the branch real on the positive axis. */
export function logGamma(z: Enclosable, opts?: number | { bits?: number }): Rect;
export function logGamma(z: CTM): CTM;
/** The Riemann zeta function. Refuses ("pole") at s = 1. */
export function zeta(s: Enclosable, opts?: number | { bits?: number }): Rect;
export function zeta(s: CTM): CTM;
/** (s − 1)ζ(s), entire. */
export function zetaTimesSMinus1(s: Enclosable, opts?: number | { bits?: number }): Rect;
export function zetaTimesSMinus1(s: CTM): CTM;
export interface ZetaZero {
  /** The interval holding the zero's imaginary part. */
  t: Bounds;
  enclosure: Rect;
  multiplicity: number;
  /** Proved to have real part exactly 1/2. */
  onCriticalLine: boolean;
  boundary: boolean;
}
/** The zeros of ζ with t1 < Im s < t2, each checked to lie on the critical line. */
export function zetaZeros(t1: RationalLike, t2: RationalLike, opts?: { bits?: number; order?: number; tolerance?: RationalLike; maxCells?: number }): { zeros: ZetaZero[]; count: number; complete: boolean; allOnCriticalLine: boolean; unresolved: Unresolved[] };

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text ("1,200", or "600, plus 800 in graphics"). */
  estimate: string;
}
/** What is still planned: nothing, every phase is built; the list is empty. */
export const plan: readonly PlannedPhase[];
