//
// Main program of the test suite for the make_ref class template.
//
//
// Original author Rob Kutschke
//

#include <iostream>

#include "Mu2eUtilities/test/makeIt.hh"

using namespace std;

int main (){

  // Const may appear in either place.
  cet::maybe_ref<mu2e::TestClass const > o1(makeItConst(1,2));
  cet::maybe_ref<const mu2e::TestClass > o2(makeItConst(3,4));
  cout << "o1: " << o1.ref() << endl;
  cout << "o2: " << o2.ref() << endl;

  // Check swap
  o2.swap(o1);
  cout << "After swap: \n"
       << "  o1: " << o1.ref()
       << "  o2: " << o2.ref()
       << endl;

  cout << "Check bool o1: " << o1 << " " << o1.isValid() << endl;

  // This should give a comiple-time error.
  //mu2e::TestClass& tc(o1.ref());

  mu2e::TestClass const& tc1(o1.ref());
  // This should give a compile-time error.
  //tc1.setDatum1(42);

  // Check the non-const maybe_ref
  cet::maybe_ref<mu2e::TestClass> o3(makeIt(5,6));
  mu2e::TestClass& tc3(o3.ref());
  cout << "tc3: Original: " << tc3 << endl;
  tc3.setDatum1(42);
  tc3.setDatum2(43);
  cout << "tc3: Revised:  " << tc3 << endl;

  // Construct an invalid object.
  cet::maybe_ref<mu2e::TestClass const> o4(makeItConst(0,6));
  cout << "Check bool o4: " << o4 << " " << o4.isValid() << endl;

  // Check that accessing an invalid object throws correctly.
  int rc = -1;
  try
    {
      cout << o4.ref() << endl;
      cout << "Failed to throw an exception!" << endl;
      rc = 1;
    }
  catch (std::logic_error& x)
    {
      cout << "Caught the expected exeption\n"
           << x.what() << endl;
      rc = 0;
    }
  catch (...)
    {
      cout << "Threw the wrong exception type! " << endl;
      rc = 2;
    }
  return rc;
}
