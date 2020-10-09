#ifndef Mu2eUtilities_test_makeIt_hh
#define Mu2eUtilities_test_makeIt_hh
//
// Part of the test suite for the make_ref class template.
//
//
// Original author Rob Kutschke
//

#include "TestTools/inc/TestClass.hh"
#include "cetlib/maybe_ref.h"

cet::maybe_ref<mu2e::TestClass const> makeItConst( int i, int j);
cet::maybe_ref<mu2e::TestClass> makeIt( int i, int j);

#endif /* Mu2eUtilities_test_makeIt_hh */
