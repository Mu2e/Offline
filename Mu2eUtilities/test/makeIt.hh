#ifndef makeIt_HH
#define makeIt_HH
//
// Part of the test suite for the make_ref class template.
//
//  $Id: makeIt.hh,v 1.1 2011/03/04 19:54:17 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/03/04 19:54:17 $
//
// Original author Rob Kutschke
//

#include "TestTools/inc/TestClass.hh"
#include "Mu2eUtilities/inc/maybe_ref.hh"

cet::maybe_ref<mu2e::TestClass const> makeItConst( int i, int j);
cet::maybe_ref<mu2e::TestClass> makeIt( int i, int j);

#endif
