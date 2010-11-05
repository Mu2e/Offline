//
// Test #1 for the MapVector class.
//
// $Id: testMapVector1.cc,v 1.1 2010/11/05 15:19:27 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/05 15:19:27 $
//
//   Original author Rob Kutschke

#include <iostream>
#include <cassert>

#include "../inc/MapVector.hh"

using namespace std;

int main(){

  // Test the default c'tor and the normal [] operator.
  MapVector<int> m1;
  m1[1] = 10;
  m1[2] = 20;

  // Test some accessors
  assert ( m1.size()         ==  2 );
  assert ( m1[1]             == 10 );
  assert ( m1.findOrThrow(2) == 20 );

  // Test the copy c'tor
  MapVector<int> m2(m1);
  assert ( m2.size()         ==  2 );
  assert ( m2[1]             == 10 );
  assert ( m2.findOrThrow(2) == 20 );

  // Modify and extend the first map.
  m1[1] = -10;
  m1[2] = -20;
  m1[3] = -30;
  assert ( m1.size()         ==   3 && "After modify" );
  assert ( m1[1]             == -10 && "After modify" );
  assert ( m1[2]             == -20 && "After modify" );
  assert ( m1[3]             == -30 && "After modify" );

  // Test the swap operator.
  m1.swap(m2);
  assert ( m1.size()         ==   2 && "After swap" );
  assert ( m2.size()         ==   3 && "After swap" );

  // Test c'tor from interators
  MapVector<int> m3(m1.begin(),m1.end());
  assert ( m3.size()         ==   2 );
  assert ( m3[1]             ==  10 );
  assert ( m3[2]             ==  20 );

  // Test a successful insertOrThrow
  MapVector<int>::iterator i = m3.insertOrThrow(make_pair<size_t,int>(4,40));
  assert ( i->first  ==  4 );
  assert ( i->second == 40 );

  // Test the const accessor
  MapVector<int> const& m4(m1);
  assert ( m4[2] ==  20 );

  // Now test things that should throw.

  // Test that insertOrThrow will throw if the key already exists.
  {
    bool pass = false;
    try {
      m1.insertOrThrow(make_pair<size_t,int>(2,20));
    } catch (std::out_of_range e ) {
      pass = true;
    }
    assert ( pass && "insertOrThrow did not throw when it should have" );
  }

  // Test that findOrThrow will throw if the key does not exist.
  {
    bool pass = false;
    try {
      double d = m1.findOrThrow(4);
    } catch (std::out_of_range e ) {
      pass = true;
    }
    assert ( pass && "findOrThrow did not throw when it should have" );
  }

  // Test that operator[](key_type key) const will throw if the key does not exist.
  {
    bool pass = false;
    try {
      double d = m4[4];
    } catch (std::out_of_range e ) {
      pass = true;
    }
    assert ( pass && "operator[](key_type) const did not throw when it should have" );
  }

  cout << "testMapVector1 passed" << endl;

  return 0;
}
