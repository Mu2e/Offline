#ifndef Sandbox_TracerProduct_hh
#define Sandbox_TracerProduct_hh
//
// A test class that makes printout whenever its methods are called.
// Each object has a "value" as a mock up of its data package plus
// a unique serial number.
//
// $Id: TracerProduct.hh,v 1.5 2013/04/07 19:41:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/04/07 19:41:12 $
//
// Original author Rob Kutschke
//

#include <ostream>

namespace mu2e {

  class TracerProduct{

  public:

    TracerProduct();
    explicit TracerProduct( int aval );
    TracerProduct( TracerProduct const& );
    TracerProduct( TracerProduct && );
    ~TracerProduct();

    TracerProduct& operator=( TracerProduct const& );
    TracerProduct& operator=( TracerProduct && );
    bool operator<( const TracerProduct& p );

    void swap( TracerProduct& );

    int val()    const { return val_; }
    int serial() const { return serial_; }

    void print ( std::ostream& ) const;

  private:

    int val_;
    int serial_;

    // Maintain a monotonically increasing serial number, starting at 0.
    // When called, return the serial number and increment it.
    static int count();

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   TracerProduct const& hit){
    hit.print(ost);
    return ost;
  }

} // namespace mu2e

#endif /* Sandbox_TracerProduct_hh */
