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

  class GammaDaughterProcessFilter : public art::EDFilter {
  public:
    explicit GammaDaughterProcessFilter(fhicl::ParameterSet const& pset);
    virtual ~GammaDaughterProcessFilter() { }

    bool filter( art::Event& event);

  private:

    // Module label of the g4 module that produced the particles
    std::string _g4ModuleLabel;
    bool _doFilter;
    int _verbosityLevel;
  };

  GammaDaughterProcessFilter::GammaDaughterProcessFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
    _doFilter(pset.get<bool>("doFilter", true)),
    _verbosityLevel(pset.get<int>("verbosityLevel", 0))
  {
  }

  bool GammaDaughterProcessFilter::filter(art::Event& event) {
    if(!_doFilter) return true;

    art::Handle<SimParticleCollection> simPCH;
    event.getByLabel(_g4ModuleLabel, simPCH);
    bool retval(true);
    const SimParticleCollection& simPC = *simPCH;
    if (_verbosityLevel >0) {
      std::cout << "SimParticleCollection has " << simPC.size() << " particles" << std::endl;
    }

    for (const auto& simPMVO : simPC) {
      const mu2e::SimParticle& simP = simPMVO.second;
      if(_verbosityLevel > 1) std::cout << "SimParticle: pdgId = " << simP.pdgId() 
					<< " end status = " << simP.endG4Status() 
					<<  " creationCode = " << simP.creationCode()
					<< " stoppingCode = " << simP.stoppingCode()
					<< std::endl;
      if(simP.stoppingCode() == ProcessCode::Mu2eGammaDaughterCut) {
	retval = false;
	break;
      }
    }
    if(_verbosityLevel > 0) {
      std::cout << "Filter returning " << retval << std::endl;
    }
    return retval;
  }
}

using mu2e::GammaDaughterProcessFilter;
DEFINE_ART_MODULE(GammaDaughterProcessFilter);
