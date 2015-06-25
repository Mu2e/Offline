#ifndef TrkPatRec_PayloadSaver_hh
#define TrkPatRec_PayloadSaver_hh
//
// Extract the persistent payload from the transient track objects.
// and put it into the event.
//
// $Id: PayloadSaver.hh,v 1.2 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Contact person Rob Kutschke
//

// Cannot use forward declaration because of typedef.
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/KalRepCollection.hh"

// Forward declarations
namespace art{
  class Event;
  class ProductID;
}

namespace fhicl{
  class ParameterSet;
}

namespace mu2e {

  class PayloadSaver{

  public:

    explicit PayloadSaver( fhicl::ParameterSet const& pset );

    void put( KalRepCollection& tracks,
              art::ProductID const& tracksID,
              art::Event&           event );

  };

} // namespace mu2e

#endif /* TrkPatRec_PayloadSaver_hh */
