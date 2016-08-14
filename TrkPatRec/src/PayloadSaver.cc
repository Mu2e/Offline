//
// Extract the persistent payload from the transient track objects.
// and put it into the event.
//
// $Id: PayloadSaver.cc,v 1.4 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Contact person Rob Kutschke
//

#include <memory>

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Provenance/ProductID.h"
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

    std::unique_ptr<KalRepPayloadCollection> payload(new KalRepPayloadCollection() );

    // Do the hard work here.
    // The ProductId is needed to make the Ptr that lives inside each KalRepPayload.

    event.put(std::move(payload));
  }

} // namespace mu2e
