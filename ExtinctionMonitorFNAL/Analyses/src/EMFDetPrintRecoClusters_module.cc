// Printout ExtMonFNAL clusers
//
// Andrei Gaponenko, 2012

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"


#include "art/Framework/Core/ModuleMacros.h"

#include "ExtinctionMonitorFNAL/Analyses/inc/GenericCollectionPrinter.hh"

namespace mu2e {

  class EMFDetPrintRecoClusters : public GenericCollectionPrinter<ExtMonFNALRecoClusterCollection> {
  public:
    explicit EMFDetPrintRecoClusters(const fhicl::ParameterSet& pset)
      : GenericCollectionPrinter(pset)
    {}

    void analyze(const art::Event& event);
  };


  //================================================================
  void EMFDetPrintRecoClusters::analyze(const art::Event& event) {
    //GenericCollectionPrinter<ExtMonFNALRecoClusterCollection>::analyze(event);
    //std::cout<<"*** Now the same info as viewed by-plane:"<<std::endl;

    std::cout<<"EMFDetPrintRecoClusters: reco pixel clusters by plane:"<<std::endl;

    art::Handle<ExtMonFNALRecoClusterCollection> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);
    const ExtMonFNALRecoClusterCollection& inputs(*ih);

    for(unsigned plane=0; plane<inputs.nplanes(); ++plane) {
      std::cout<<"   Plane "<<plane<<std::endl;
      ExtMonFNALRecoClusterCollection::PlaneClusters pc = inputs.clusters(plane);
      for(unsigned i=0; i<pc.size(); ++i) {
        std::cout<<"       "<<pc[i]<<std::endl;
      }
    }
  }

}

DEFINE_ART_MODULE(mu2e::EMFDetPrintRecoClusters);
