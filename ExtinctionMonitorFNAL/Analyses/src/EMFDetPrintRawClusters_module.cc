// Printout ExtMonFNAL raw clusters
//
// Andrei Gaponenko, 2012

#include "Offline/RecoDataProducts/inc/ExtMonFNALRawCluster.hh"


#include "Offline/ExtinctionMonitorFNAL/Analyses/inc/GenericCollectionPrinter.hh"

namespace mu2e {

  class EMFDetPrintRawClusters : public GenericCollectionPrinter<ExtMonFNALRawClusterCollection> {
  public:
    explicit EMFDetPrintRawClusters(const fhicl::ParameterSet& pset)
      : GenericCollectionPrinter(pset)
    {}

  };

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetPrintRawClusters)
