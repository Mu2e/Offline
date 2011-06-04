#ifndef Sandbox_TransientProduct00_hh
#define Sandbox_TransientProduct00_hh
//
// A test data class that contains a bare pointer.
// Used for tests of making transient-only data products.
//
// $Id: TransientProduct00.hh,v 1.1 2011/06/04 18:00:36 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/04 18:00:36 $
//
// Original author Rob Kutschke
//

#include "RecoDataProducts/inc/StrawHit.hh"

namespace mu2e {

  //  class StrawHit;

  class TransientProduct00{

  public:

    TransientProduct00():
      hit_(0){
    }

    TransientProduct00( StrawHit const& hit):
      hit_(&hit){
    }
    // Accept compiler generated d'tor, copy c'tor and operator=.

    // Accessors
    StrawHit const& hit() const { return *hit_;}

  private:

    StrawHit const* hit_;

  };

} // namespace mu2e

#endif /* Sandbox_TransientProduct00_hh */
