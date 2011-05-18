//
// A trivial class instrumented with some printout.  To
// be used in tests of other classes.
//
//  $Id: TestClass.cc,v 1.2 2011/05/18 02:27:19 wb Exp $
//  $Author: wb $
//  $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include <iostream>
#include <iomanip>

#include "TestTools/inc/TestClass.hh"

using namespace std;

namespace mu2e {

  int TestClass::counter = -1;

  TestClass::TestClass ( ):
    _datum1(0),
    _datum2(0){
    _serial = ++counter;
    cout << "Default constructor: "
         << setw(2)
         << _datum1 << " "
         << _datum2
         << "  _serial: "
         << _serial
         << endl;
  }

  TestClass::TestClass ( int i, int j ):
    _datum1(i),
    _datum2(j){
    _serial = ++counter;
    cout << "Int constructor:     "
         << setw(2)
         << _datum1 << " "
         << _datum2
         << "  _serial: "
         << _serial
         << endl;
  }

  TestClass::TestClass ( const TestClass& m ):
    _datum1(m._datum1),
    _datum2(m._datum2){
    _serial = ++counter;
    cout << "Copy constructor:    "
         << setw(2)
         << _datum1 << " "
         << _datum2
         << "  From _serial: "
         << m._serial
         << " To _serial: "
         << _serial
         << endl;
  }

  TestClass::~TestClass ( ){
    cout << "Destructor:          "
         << setw(2)
         << _datum1 << " "
         << _datum2
         << "  _serial: "
         << _serial
         << endl;
  }

  TestClass& TestClass::operator=(TestClass rhs){
    cout << "Invoking assignment:  From: "
         << rhs._serial
         << "  To: "
         << _serial
         << endl;
    rhs.swap(*this);
    return *this;
  }

} // end namespace mu2e
