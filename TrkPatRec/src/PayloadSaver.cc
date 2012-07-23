//
// Extract the persistent payload from the transient track objects.
// and put it into the event.
//
// $Id: PayloadSaver.cc,v 1.2 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Contact person Rob Kutschke
//

#include <memory>

#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Provenance/ProductID.h"
#include "fhiclcpp/ParameterSet.h"

#include "TrkPatRec/inc/PayloadSaver.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"

namespace mu2e {

  PayloadSaver::PayloadSaver( fhicl::ParameterSet const& pset ){
    // Someday we will have run-time config information.  Extract that info here.
  }

  void PayloadSaver::put( KalRepCollection& tracks,
                          art::ProductID const& tracksID,
                          art::Event& event ){

    std::auto_ptr<KalRepPayloadCollection> payload(new KalRepPayloadCollection() );

    // Do the hard work here.
    // The ProductId is needed to make the Ptr that lives inside each KalRepPayload.

    event.put(payload);
  }

} // namespace mu2e
