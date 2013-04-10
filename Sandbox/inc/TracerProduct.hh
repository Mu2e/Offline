#ifndef Sandbox_TracerProduct_hh
#define Sandbox_TracerProduct_hh
//
// A test class that makes printout whenever its methods are called.
// Each object has a "value" as a mock up of its data package plus
// a unique serial number.
//
// $Id: TracerProduct.hh,v 1.6 2013/04/10 14:38:02 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/04/10 14:38:02 $
//
// Original author Rob Kutschke
//

#include <ostream>

namespace mu2e {

  class TracerProduct{

  public:

    TracerProduct();
    ~TracerProduct();

#ifndef __GCCXML__
    explicit TracerProduct( int aval );
    TracerProduct( TracerProduct const& );
    TracerProduct( TracerProduct && );

    TracerProduct& operator=( TracerProduct const& );
    TracerProduct& operator=( TracerProduct && );
    bool operator<( const TracerProduct& p );

    void swap( TracerProduct& );

    int val()    const { return val_; }
    int serial() const { return serial_; }

    void print ( std::ostream& ) const;
#endif

  private:

    int val_;
    int serial_;

#ifndef __GCCXML__
    // Maintain a monotonically increasing serial number, starting at 0.
    // When called, return the serial number and increment it.
    static int count();
#endif

  };

#ifndef __GCCXML__
  inline std::ostream& operator<<( std::ostream& ost,
                                   TracerProduct const& hit){
    hit.print(ost);
    return ost;
  }
#endif

} // namespace mu2e

#endif /* Sandbox_TracerProduct_hh */
