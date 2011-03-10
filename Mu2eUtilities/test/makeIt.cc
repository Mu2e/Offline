//
// Part of the test suite for the make_ref class template.
//
//  $Id: makeIt.cc,v 1.2 2011/03/10 00:00:58 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/03/10 00:00:58 $
//
// Original author Rob Kutschke
//

#include <list>
#include "Mu2eUtilities/test/makeIt.hh"


// A cache to make sure all of the test objects live until
// the end of the job.
static std::list<mu2e::TestClass> v;

cet::maybe_ref<mu2e::TestClass> makeIt( int i, int j){
  if ( i == 0 ){
    return cet::maybe_ref<mu2e::TestClass>();
  }
  v.push_back(mu2e::TestClass(i,j));
  return cet::maybe_ref<mu2e::TestClass>(v.back());
}

cet::maybe_ref<mu2e::TestClass const > makeItConst( int i, int j){
  if ( i == 0 ){
    return cet::maybe_ref<mu2e::TestClass const >();
  }
  v.push_back(mu2e::TestClass(i,j));
  return cet::maybe_ref<mu2e::TestClass const >(v.back());
}
