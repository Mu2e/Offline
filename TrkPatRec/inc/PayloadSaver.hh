#ifndef TrkPatRec_PayloadSaver_hh
#define TrkPatRec_PayloadSaver_hh
//
// Extract the persistent payload from the transient track objects.
// and put it into the event.
//
// $Id: PayloadSaver.hh,v 1.1 2012/07/03 04:19:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 04:19:14 $
//
// Contact person Rob Kutschke
//

// Cannot use forward declaration because of typedef.
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh"

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

    void put( TrkRecoTrkCollection& tracks,
              art::ProductID const& tracksID,
              art::Event&           event );

  };

} // namespace mu2e

#endif /* TrkPatRec_PayloadSaver_hh */
