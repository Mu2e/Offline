#ifndef Mu2eUtilities_test_makeIt_hh
#define Mu2eUtilities_test_makeIt_hh
//
// Part of the test suite for the make_ref class template.
//
//  $Id: makeIt.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
//  $Author: greenc $
//  $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include "TestTools/inc/TestClass.hh"
#include "Mu2eUtilities/inc/maybe_ref.hh"

cet::maybe_ref<mu2e::TestClass const> makeItConst( int i, int j);
cet::maybe_ref<mu2e::TestClass> makeIt( int i, int j);

#endif /* Mu2eUtilities_test_makeIt_hh */
