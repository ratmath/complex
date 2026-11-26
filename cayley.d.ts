// Type declarations for @ratmath/complex/cayley.

import type { Rational } from "@ratmath/core";
import type { Field, RationalLike } from "./index.js";

/** The rationals as a field object. */
export const QQ: Field<Rational>;

export interface AlgebraProperties<T = Rational> {
  name: string;
  dimension: number;
  commutative: boolean;
  associative: boolean;
  alternative: boolean;
  flexible: boolean;
  /** N(xy) = N(x)N(y). */
  composition: boolean;
  /** null, or two nonzero elements whose product is 0. */
  zeroDivisors: [Hyper<T>, Hyper<T>] | null;
}

/** One level of the doubling construction: an algebra of dimension 2^level over a field. */
export class CayleyDickson<T = Rational> {
  readonly level: number;
  readonly dim: number;
  readonly field: Field<T>;
  readonly name: string;
  readonly labels: readonly string[];
  readonly zero: Hyper<T>;
  readonly one: Hyper<T>;
  readonly units: Hyper<T>[];
  constructor(level: number, field?: Field<T>);
  basis(i: number): Hyper<T>;
  of(...coeffs: unknown[]): Hyper<T>;
  parse(text: string): Hyper<T>;
  properties(): AlgebraProperties<T>;
  toString(): string;
}

/** Anything an element's methods accept as the other operand. */
export type HyperLike<T = Rational> = Hyper<T> | RationalLike | readonly unknown[];

export class Hyper<T = Rational> {
  readonly algebra: CayleyDickson<T>;
  readonly coeffs: readonly T[];
  readonly real: T;
  readonly imaginary: Hyper<T>;
  add(o: HyperLike<T>): Hyper<T>;
  sub(o: HyperLike<T>): Hyper<T>;
  mul(o: HyperLike<T>): Hyper<T>;
  neg(): Hyper<T>;
  scale(c: unknown): Hyper<T>;
  conj(): Hyper<T>;
  norm(): T;
  inv(): Hyper<T>;
  /** The x with x·b = this. */
  rightDivide(b: HyperLike<T>): Hyper<T>;
  /** The x with b·x = this. */
  leftDivide(b: HyperLike<T>): Hyper<T>;
  pow(n: number): Hyper<T>;
  commutator(o: HyperLike<T>): Hyper<T>;
  associator(b: HyperLike<T>, c: HyperLike<T>): Hyper<T>;
  isZero(): boolean;
  equals(o: HyperLike<T>): boolean;
  toString(): string;
  toJSON(): string[];
}

export function cayleyDickson(level: number): CayleyDickson<Rational>;
export function cayleyDickson<T>(level: number, field: Field<T>): CayleyDickson<T>;
export const quaternions: CayleyDickson<Rational>;
export const octonions: CayleyDickson<Rational>;
export const sedenions: CayleyDickson<Rational>;
export function quaternion(w?: RationalLike, x?: RationalLike, y?: RationalLike, z?: RationalLike): Hyper;
export function octonion(...coeffs: RationalLike[]): Hyper;

/** A quaternion, its four coordinates, or text. */
export type QuaternionLike = Hyper | RationalLike | readonly RationalLike[];
export type Vector3 = readonly [RationalLike, RationalLike, RationalLike];

export function rotationMatrix(q: QuaternionLike): Rational[][];
export function rotate(q: QuaternionLike, v: Vector3): Rational[];
export function fromAxis(axis: Vector3, t?: RationalLike): Hyper;
export function rotationOf(q: QuaternionLike): { axis: Rational[]; cos: Rational; tanHalfSquared: Rational | null };
export function primitive(q: QuaternionLike): Hyper;
export function fromRotationMatrix(M: readonly (readonly RationalLike[])[]): Hyper;
export function interpolate(q0: QuaternionLike, q1: QuaternionLike, t: RationalLike): Hyper;

export type HurwitzLike = Hurwitz | Hyper | readonly RationalLike[];
/** A quaternion whose coordinates are all whole numbers or all halves of odd numbers. */
export class Hurwitz {
  /** Doubled coordinates: [1n, 1n, 1n, 1n] is (1 + i + j + k)/2. */
  readonly d: readonly bigint[];
  constructor(doubled: readonly (bigint | number)[]);
  static from(x: HurwitzLike | RationalLike, ...rest: RationalLike[]): Hurwitz;
  readonly isLipschitz: boolean;
  norm(): bigint;
  conj(): Hurwitz;
  neg(): Hurwitz;
  add(o: HurwitzLike): Hurwitz;
  sub(o: HurwitzLike): Hurwitz;
  mul(o: HurwitzLike): Hurwitz;
  /** this = quotient·b + remainder with N(remainder) < N(b). */
  divmod(b: HurwitzLike): { quotient: Hurwitz; remainder: Hurwitz };
  isZero(): boolean;
  equals(o: HurwitzLike): boolean;
  toQuaternion(): Hyper;
  toIntegers(): bigint[];
  toString(): string;
}
export function hurwitz(...coordinates: (RationalLike | HurwitzLike)[]): Hurwitz;
export function hurwitzUnits(): Hurwitz[];
export function hurwitzGcd(a: HurwitzLike, b: HurwitzLike): Hurwitz;
export function quaternionOfNorm(p: bigint | number): Hurwitz;
export function fourSquares(n: bigint | number): bigint[];
export function fourSquaresCount(n: bigint | number): bigint;
