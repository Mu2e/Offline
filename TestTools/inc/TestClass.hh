#ifndef TESTCLASS_HH
#define TESTCLASS_HH
//
// A trivial class instrumented with some printout.  To
// be used in tests of other classes.
//
//  $Id: TestClass.hh,v 1.1 2011/03/04 19:54:17 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/03/04 19:54:17 $
//
// Original author Rob Kutschke
//

#include <algorithm>
#include <ostream>

namespace mu2e {

  class TestClass{

  public:
    TestClass ( );
    TestClass ( int i, int j=-1 );
    TestClass ( const TestClass& m );
    ~TestClass();

    void swap ( TestClass& a ){
      std::swap( *this, a );
    }

    TestClass& operator=(TestClass rhs);

    int datum1() const { return _datum1;}
    int datum2() const { return _datum2;}
    int serial() const { return _serial;}

    void setDatum1( int i) {
      _datum1 = i;
    }

    void setDatum2( int j) {
      _datum2 = j;
    }

  private:

    int _datum1;
    int _datum2;

    // Give a unique serial number to each object.
    static int counter;
    int _serial;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   TestClass const& a){
    ost << "serial:   "  << a.serial()
        << "  _datum1: " << a.datum1()
        << "  _datum2: " << a.datum2();
    return ost;
  }

} // namespace mu2e

#endif
