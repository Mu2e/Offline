
// ======================================================================
//
// GenParticleMomFilter_module: allows filtering on the energy
//   of the GenParticle
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"



#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class GenParticleMomFilter : public art::EDFilter {
    public:
      explicit GenParticleMomFilter(fhicl::ParameterSet const& pset);

    private:
      bool beginRun(art::Run& run) override;
      bool endRun(art::Run& run) override;
      bool filter(art::Event& event) override;
      
      art::InputTag _genParticleModule;
      double _momCutoff;
      PDGCode::type _cutoffPDG;
      GenId _cutoffGenId;
      unsigned      _nevt, _npass;
  };

  GenParticleMomFilter::GenParticleMomFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _genParticleModule(pset.get<std::string>("genParticleModule","compressDigiMCs")),
    _momCutoff(pset.get<double>("MomentumCutoff")),
    _cutoffPDG(PDGCode::type(pset.get<int>("CutoffPDG",0))),
    _cutoffGenId(GenId::findByName(pset.get<std::string>("CutoffGenId",""))),
    _nevt(0), _npass(0){}

  bool GenParticleMomFilter::beginRun(art::Run& run) {
    return true;
  }
  bool GenParticleMomFilter::endRun(art::Run& run) {
    return true;
  }

  bool GenParticleMomFilter::filter(art::Event& event) {
    ++_nevt;

    auto genColl = event.getValidHandle<GenParticleCollection>( _genParticleModule);

    // find highest momentum gen particle that passes cuts
    double mom = 0;
    for ( const auto& i: *genColl ) {
      if ((i.pdgId() == _cutoffPDG || _cutoffPDG == PDGCode::null) && i.generatorId() == _cutoffGenId) {
        if (i.momentum().vect().mag() > mom)
          mom = i.momentum().vect().mag();
      }
    }
    if (mom < _momCutoff){
      return false;
    }

    ++_npass;
    return true;
  }
}

using mu2e::GenParticleMomFilter;
DEFINE_ART_MODULE(GenParticleMomFilter);
