//
// Part of the test suite for the make_ref class template.
//
//  $Id: makeIt.cc,v 1.1 2011/03/04 19:54:17 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/03/04 19:54:17 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/test/makeIt.hh"

// Yes, this leaks.  Not important for this test.
cet::maybe_ref<mu2e::TestClass> makeIt( int i, int j){
  if ( i == 0 ){
    return cet::maybe_ref<mu2e::TestClass>();
  }
  return cet::maybe_ref<mu2e::TestClass>(new mu2e::TestClass(i,j));
}

cet::maybe_ref<mu2e::TestClass const > makeItConst( int i, int j){
  if ( i == 0 ){
    return cet::maybe_ref<mu2e::TestClass const >();
  }
  return cet::maybe_ref<mu2e::TestClass const >(new mu2e::TestClass(i,j));
}
