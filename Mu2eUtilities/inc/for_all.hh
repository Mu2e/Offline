#ifndef Mu2eUtilities_for_all_hh
#define Mu2eUtilities_for_all_hh

//
// A free function to iterate over all elements of a container
// and call a function, with two arguments:
//  1) A reference to the element in the container ( non-const)
//  2) An arbitrary second object.
// A later version of C++ will allow a ... notation for an arbitrary
// number of additional objects.
//
// $Id: for_all.hh,v 1.4 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//

namespace mu2e {

  template < class C, class F, class Acc>
  inline void for_all ( C const& c, F f, Acc& acc){
    for ( typename C::const_iterator i=c.begin(), e=c.end();
          i != e; ++i ){
      f(*i,acc);
    }
  }

  // A related version in which it is possible to modify the members of
  // the collection.
  template < class C, class F, class Acc>
  inline void for_all ( C& c, F f, Acc& acc){
    for ( typename C::const_iterator i=c.begin(), e=c.end();
          i != e; ++i ){
      f(*i,acc);
    }
  }

}

#endif /* Mu2eUtilities_for_all_hh */
