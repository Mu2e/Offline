//
//

// Mu2e includes.
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class GammaDaughterFilter : public art::EDFilter {
  public:
    explicit GammaDaughterFilter(fhicl::ParameterSet const& pset);
    virtual ~GammaDaughterFilter() { }

    bool filter( art::Event& event);

  private:

    // Module label of the g4 module that produced the particles
    std::string _g4ModuleLabel;
    bool _doFilter;
    int _verbosityLevel;
  };

  GammaDaughterFilter::GammaDaughterFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
    _doFilter(pset.get<bool>("doFilter", true)),
    _verbosityLevel(pset.get<int>("verbosityLevel", 0))
  {
  }

  bool GammaDaughterFilter::filter(art::Event& event) {
    if(!_doFilter) return true;

    art::Handle<SimParticleCollection> simPCH;
    event.getByLabel(_g4ModuleLabel, simPCH);
    bool retval(true);
    const SimParticleCollection& simPC = *simPCH;
    if (_verbosityLevel >0) {
      cout << "SimParticleCollection has " << simPC.size() << " particles" << endl;
    }

    for (const auto& simPMVO : simPC) {
      const mu2e::SimParticle& simP = simPMVO.second;
      if(_verbosityLevel > 1) cout << "SimParticle: pdgId = " << simP.pdgId() 
				   << " end status = " << simP.endG4Status() 
				   <<  " creationCode = " << simP.creationCode().id()
				   << " stoppingCode = " << simP.stoppingCode().id()
				   << endl;
      if(simP.stoppingCode() == ProcessCode::mu2eLowEnergyGammaKilled) {
	retval = false;
	break;
      }
    }
    return retval;
  }
}

using mu2e::GammaDaughterFilter;
DEFINE_ART_MODULE(GammaDaughterFilter);
