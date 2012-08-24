// Printout ExtMonFNAL raw hits and associated truth
//
// Andrei Gaponenko, 2012

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  //================================================================
  class EMFDetPrintTruth : public art::EDAnalyzer {
    std::string _inModuleLabel;
    std::string _inInstanceName;
  public:
    explicit EMFDetPrintTruth(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetPrintTruth::EMFDetPrintTruth(const fhicl::ParameterSet& pset)
    : _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
    , _inInstanceName(pset.get<std::string>("inputInstanceName"))
  {}

  //================================================================
  void EMFDetPrintTruth::analyze(const art::Event& event) {

    art::Handle<ExtMonFNALRawHitCollection> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);

    const ExtMonFNALRawHitCollection& inputs(*ih);

    art::FindManyP<SimParticle,ExtMonFNALHitTruthBits> r2t(ih, event, _inModuleLabel);

    std::cout<<"EMFDetPrintTruth: inModuleLabel = "<<_inModuleLabel<<", inInstanceName = "<<_inInstanceName<<std::endl;

    for(ExtMonFNALRawHitCollection::const_iterator i=inputs.begin(); i!=inputs.end(); ++i) {
      std::cout<<"event "<<event.id()<<", hit "<<*i<<std::endl;

      std::vector<art::Ptr<SimParticle> > particles;
      std::vector<const ExtMonFNALHitTruthBits*> charges;

      r2t.get((i-inputs.begin()), particles, charges);
      std::cout<<"got truth: num particles = "<<particles.size()<<", num charges "<<charges.size()<<std::endl;;
      for(unsigned ip=0; ip<particles.size(); ++ip) {
        std::cout<<"track "<<particles[ip]->id()<<", charge "<<charges[ip]->charge()<<std::endl;
      }
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetPrintTruth);
