// Printout ExtMonFNAL hits
//
// Andrei Gaponenko, 2012

#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"


#include "Offline/ExtinctionMonitorFNAL/Analyses/inc/GenericCollectionPrinter.hh"

namespace mu2e {

  class EMFDetPrintSim : public GenericCollectionPrinter<ExtMonFNALSimHitCollection> {
  public:
    explicit EMFDetPrintSim(const fhicl::ParameterSet& pset)
      : GenericCollectionPrinter(pset)
    {}

  };

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetPrintSim)
