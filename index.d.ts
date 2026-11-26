// Type declarations for @ratmath/complex.

import type { Rational, RationalInterval } from "@ratmath/core";
import type { Refinement } from "@ratmath/core/refine";

/** Anything read as a real rational: a Rational, bigint, number, or text in Core's notation. */
export type RationalLike = Rational | bigint | number | string;
/** Anything read as an exact complex number. */
export type ComplexLike = Complex | RationalLike | readonly [RationalLike, RationalLike] | { re: RationalLike; im: RationalLike };
/** An exact number or an enclosure. */
export type Enclosable = ComplexLike | Rect | Disc;

export type ComplexErrorCode = "division-by-zero" | "branch-cut" | "branch-point" | "pole" | "precision" | "bad-option" | "no-quotient" | "not-a-rotation" | "too-large";
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
export function sin(z: Enclosable, opts?: Precision): Rect;
export function cos(z: Enclosable, opts?: Precision): Rect;
export function tan(z: Enclosable, opts?: Precision): Rect;
export function sinh(z: Enclosable, opts?: Precision): Rect;
export function cosh(z: Enclosable, opts?: Precision): Rect;
export function tanh(z: Enclosable, opts?: Precision): Rect;

export function log(z: Enclosable, opts?: number | BranchOptions): Branched;
export function sqrt(z: Enclosable, opts?: number | BranchOptions): Branched;
export function root(z: Enclosable, n: number, opts?: number | BranchOptions): Branched;
export function pow(z: Enclosable, w: Enclosable, opts?: number | BranchOptions): Branched;
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

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text ("1,200", or "600, plus 800 in graphics"). */
  estimate: string;
}
/** What is still planned: phases 3 to 6. */
export const plan: readonly PlannedPhase[];
