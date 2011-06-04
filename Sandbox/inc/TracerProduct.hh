#ifndef Sandbox_TracerProduct_hh
#define Sandbox_TracerProduct_hh
//
// A test class that makes printout whenever its methods are called.
//
// $Id: TracerProduct.hh,v 1.1 2011/06/04 20:36:59 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/04 20:36:59 $
//
// Original author Rob Kutschke
//

namespace mu2e {

  class TracerProduct{

  public:

    TracerProduct();
    explicit TracerProduct( int aval );
    explicit TracerProduct( TracerProduct const& );
    ~TracerProduct();
    
    TracerProduct& operator=( TracerProduct const& );
    bool operator<( const TracerProduct& p );

    int val()    const { return val_; }
    int serial() const { return serial_; }
    
  private:

    int val_;
    int serial_;
    
    static int count();

  };

} // namespace mu2e

#endif /* Sandbox_TracerProduct_hh */
