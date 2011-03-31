// Angle classes, with automatic range checking.
// SignedAngle ranges from -PI to PI, UnsignedAngle from 0 to 2PI.
// The two types of angle convert implicitly into double.
//
// The Angle classes reimplement (inline) most arithmetic operators,
// instead of relying on the implicit conversion to double and converting
// the result back when it is assigned to an Angle variable. There are two
// reasons for this:
// 1) Temporary values are kept in range (example: std::cout << alpha+beta)
// 2) Optimization (the difference of two Angles is guaranteed to lie less
//    than 2PI outside the allowed range, for example)
//
// The only exception is multiplication: the "2" in "2 * alpha" could be
// interpreted as a length. If you want the result of a multiplication to be
// bounded, assign it to a temporary Angle, or use *=.
//
// Due to operator ambiguity, some mixed-type arithmetics might fail:
// for example, (int + Angle) or (Angle / long)
// The most common operations are supported:
// double +- angle; angle +-/ double; angle / int
//
// Comparison between Angles is done via conversion to double. Make sure
// you do not compare a SignedAngle to an UnsignedAngle.
//
// Since the bounds and return types for the two classes are different,
// operators cannot be factored in a single superclass.
//
// In C++, implicit conversion operators are shallow - they are not nested.
// Therefore, if you want to convert a SignedAngle into an UnsignedAngle,
// you need an intermediate cast into double. For example:
//     SignedAngle s; UnsignedAngle u = double(s);
// Adding a direct conversion between the two types of angle would result
// in several ambiguities.
//
// Author: Paolo Gatti, University of Padova / INFN

#ifndef ANGLE_HH
#define ANGLE_HH

#include <iostream>
#include <cmath>

#ifdef __GNUG__
#define inline __attribute((always_inline))
#endif

//********************************************************************
// Signed angles, ranging from -PI to PI
class SignedAngle {
  double _value;

public:
// Range check/fix functions
  static inline double fixRange(const double v) {
    double shift = floor(v/(2*M_PI)+0.5);
    return (v - 2*M_PI*shift);
  }
// Constructor
  inline SignedAngle(double angle = 0.0) {
    _value = fixRange(angle);
  }
// No destructor; use standard copy constructor and assignment operator

// Implicit conversion to double
  inline operator double() const {return _value;}

// Self-modifying operations
  inline SignedAngle& operator += (const SignedAngle& ang) {
    _value = fixRange(_value+ang._value);
    return *this;
  }
  inline SignedAngle& operator -= (const SignedAngle& ang) {
    _value = fixRange(_value-ang._value);
    return *this;
  }
  inline SignedAngle& operator *= (double dbl) {
    _value = fixRange(_value*dbl);
    return *this;
  }
  inline SignedAngle& operator /= (double dbl) {
    _value = fixRange(_value/dbl);
    return *this;
  }

// Other operations
  inline friend SignedAngle operator + (double dbl, const SignedAngle& ang) {
    return SignedAngle(dbl+ang._value);
  }
  inline friend SignedAngle operator - (double dbl, const SignedAngle& ang) {
    return SignedAngle(dbl-ang._value);
  }
  inline SignedAngle operator + (double dbl) const {
    return SignedAngle(_value+dbl);
  }
  inline SignedAngle operator - (double dbl) const {
    return SignedAngle(_value-dbl);
  }
  inline SignedAngle operator + (SignedAngle ang) const {
    return SignedAngle(_value+ang._value);
  }
  inline SignedAngle operator - (SignedAngle ang) const {
    return SignedAngle(_value-ang._value);
  }
  inline SignedAngle operator / (double dbl) const {
    return SignedAngle(_value/dbl);
  }
  inline SignedAngle operator / (int dbl) const  {
    return SignedAngle(_value/dbl);
  }

  // I/O
  friend std::ostream& operator << (std::ostream& out, const SignedAngle& me) {
    out << me._value;
    return out;
  }

  friend std::istream& operator >> (std::istream& in, SignedAngle& me) {
    double v;
    in >> v;
    me._value = fixRange(v);
    return in;
  }
};

//********************************************************************
// Unsigned angles, ranging from 0 to 2PI
class UnsignedAngle {
  double _value;

public:
// Range check/fix functions
  static inline double fixRange(const double v) {
    double shift = floor(v/(2*M_PI));
    return (v - 2*M_PI*shift);
  }
// Constructor
  inline UnsignedAngle(double angle = 0.0) {
    _value = fixRange(angle);
  }
// No destructor; use standard copy constructor and assignment operator

// Implicit conversion to double
  inline operator double() const {return _value;}

// Self-modifying operations
  inline UnsignedAngle& operator += (const UnsignedAngle& ang) {
    _value = fixRange(_value+ang._value);
    return *this;
  }
  inline UnsignedAngle& operator -= (const UnsignedAngle& ang) {
    _value = fixRange(_value-ang._value);
    return *this;
  }
  inline UnsignedAngle& operator *= (double dbl) {
    _value = fixRange(_value*dbl);
    return *this;
  }
  inline UnsignedAngle& operator /= (double dbl) {
    _value = fixRange(_value/dbl);
    return *this;
  }

// Other operations
  inline friend UnsignedAngle operator + (double dbl, const UnsignedAngle& ang) {
    return UnsignedAngle(dbl+ang._value);
  }
  inline friend UnsignedAngle operator - (double dbl, const UnsignedAngle& ang) {
    return UnsignedAngle(dbl-ang._value);
  }
  inline UnsignedAngle operator + (double dbl) const {
    return UnsignedAngle(_value+dbl);
  }
  inline UnsignedAngle operator - (double dbl) const {
    return UnsignedAngle(_value-dbl);
  }
  inline UnsignedAngle operator + (UnsignedAngle ang) const {
    return UnsignedAngle(_value+ang._value);
  }
  inline UnsignedAngle operator - (UnsignedAngle ang) const {
    return UnsignedAngle(_value-ang._value);
  }
  inline UnsignedAngle operator / (double dbl) const {
    return UnsignedAngle(_value/dbl);
  }
  inline UnsignedAngle operator / (int dbl) const  {
    return UnsignedAngle(_value/dbl);
  }

  // I/O
  friend std::ostream& operator << (std::ostream& out, const UnsignedAngle& me) {
    out << me._value;
    return out;
  }

  friend std::istream& operator >> (std::istream& in, UnsignedAngle& me) {
    double v;
    in >> v;
    me._value = fixRange(v);
    return in;
  }
};

#ifdef __GNUG__
#undef inline
#endif

// By default, "Angles" are unsigned (CDF convention: 0-2PI range).
typedef UnsignedAngle Angle;

#endif /* ANGLE_HH */
