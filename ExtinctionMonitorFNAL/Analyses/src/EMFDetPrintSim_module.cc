// Printout ExtMonFNAL hits
//
// Andrei Gaponenko, 2012

#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

#include "art/Framework/Core/ModuleMacros.h"

#include "ExtinctionMonitorFNAL/Analyses/inc/GenericCollectionPrinter.hh"

namespace mu2e {

  class EMFDetPrintSim : public GenericCollectionPrinter<ExtMonFNALSimHitCollection> {
  public:
    explicit EMFDetPrintSim(const fhicl::ParameterSet& pset)
      : GenericCollectionPrinter(pset)
    {}

  };

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetPrintSim);
