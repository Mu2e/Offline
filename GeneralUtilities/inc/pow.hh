// ======================================================================
//
// Minimum-multiplication algorithm to compute powers known at compile time.
//
// Author:  Walter E Brown <wb@fnal.gov>
// Date:    2009-10-05
//
// Usage:
//   For any type T for which multiply operator is defined:
//     unsigned long const n;  // Must be initialized at compile time.
//     T a;                    // May be determined as late as run time.
//     T answer = pow<n>(a);
//     T answer = square(a);
//     T answer = cube(a);
//     T answer = fourth(a);
//
// $Id: pow.hh,v 1.1 2009/10/06 23:49:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/06 23:49:30 $
//
// ======================================================================


#ifndef COMPILE_TIME_POW
#define COMPILE_TIME_POW


// ----------------------------------------------------------------------
// Public contents declaration
// ----------------------------------------------------------------------

template< unsigned long  N
        , class          T
        >
T  pow( T x );


template< class T >
T  square( T x );

template< class T >
T  cube( T x );

template< class T >
T  fourth( T x );


// ----------------------------------------------------------------------
// Private
// ----------------------------------------------------------------------


namespace pow_details{

  // --------------------------------------------------------------------

  template< unsigned long  N >
  struct is_odd  { static  bool const  value = (N & 1UL) == 1UL; };

  // --------------------------------------------------------------------

  template< bool, class T, class >
  struct conditional             { typedef  T  type; };

  template< class T, class F >
  struct conditional<false,T,F>  { typedef  F  type; };

  // --------------------------------------------------------------------

  template< unsigned long  N
          , class          T
          >
  struct Power;

  // --------------------------------------------------------------------

  template< unsigned long  N
          , class          T
          >
  struct EvenPower
  {
    T  operator () ( T x )  { return Power< N/2UL, T >()( x * x ); }
  };

  template< class T >
  struct EvenPower<2UL,T>
  {
    T  operator () ( T x )  { return x * x; }
  };

  template< class T >
  struct EvenPower<0UL,T>
  {
    T  operator () ( T x )  { return static_cast<T>(1); }
  };

  // --------------------------------------------------------------------

  template< unsigned long  N
          , class          T
          >
  struct OddPower
  {
    T  operator () ( T x )  { return x * EvenPower< N-1UL, T >()( x ); }
  };

  template< class T >
  struct OddPower<1UL,T>
  {
    T  operator () ( T x )  { return x; }
  };

  // --------------------------------------------------------------------

  template< unsigned long  N
          , class          T
          >
  struct Power
    : public conditional< is_odd< N >::value
                        , OddPower < N, T >
                        , EvenPower< N, T >
                        >::type
  { };

  // --------------------------------------------------------------------

}  // namespace pow_details


// ----------------------------------------------------------------------
// Public definitions
// ----------------------------------------------------------------------

template< unsigned long N
        , class         T  // must provide operator *
        >
inline
T  pow( T x )  { return pow_details::Power<N,T>()( x ); }


template< class T >
inline
T  square( T x )  { return pow<2UL>( x ); }

template< class T >
inline
T  cube( T x )  { return pow<3UL>( x ); }

template< class T >
inline
T  fourth( T x )  { return pow<4UL>( x ); }


// ======================================================================
//

#endif  // COMPILE_TIME_POW

