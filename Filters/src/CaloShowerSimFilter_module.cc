//
// Select events with a minimum energy sum from CaloShowerSimfrom a single
// particle, and save those particles as 'primaries'.  Particle can also be
// selected based on PDG code
//  original author: David Brown (LBNL)
//

// Mu2e includes.
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
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

  class CaloShowerSimFilter : public art::EDFilter {
  public:
    explicit CaloShowerSimFilter(fhicl::ParameterSet const& pset);
    virtual ~CaloShowerSimFilter() { }

    bool filter( art::Event& event);

  private:

    unsigned minndigi_;
    double minpe_;
    double minetot_;
    std::vector<PDGCode::type> pdgs_;
    int diag_, debug_;
    art::InputTag cssTag_;

  };

  CaloShowerSimFilter::CaloShowerSimFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    minpe_(pset.get<double>("MinParticleEnergy")),
    minetot_(pset.get<double>("MinTotalEnergy")),
    diag_(pset.get<int>("diagLevel",0)),
    debug_(pset.get<int>("debugLevel",0)),
    cssTag_(pset.get<art::InputTag>("CaloShowerSimCollection","CaloShowerStepROFromShowerStep"))
  {
    produces<SimParticlePtrCollection>();
    auto pdgs = pset.get<std::vector<int>>("particleTypes");
    for(auto pdg : pdgs )
      pdgs_.push_back(PDGCode::type(pdg));
  }

  bool CaloShowerSimFilter::filter(art::Event& evt) {
    std::unique_ptr<SimParticlePtrCollection> output(new SimParticlePtrCollection());
    bool retval(false);
    auto cssH = evt.getValidHandle<CaloShowerSimCollection>(cssTag_);
    const CaloShowerSimCollection *csscol = cssH.product();
    // keep count of energy deposited by specific particle
    std::map<art::Ptr<SimParticle>,float> emap;
    double etot(0.0);
    for(auto const& css : *csscol) {
      art::Ptr<SimParticle> const& sp = css.sim();
      double energy = css.energy();
      etot += energy;
      if(debug_ > 0)std::cout <<"SimParticle PDG = " << sp->pdgId() 
      << " Crystal " << css.crystalId()
      << " Shower Energy = " << energy << std::endl;
      auto pdgfnd = std::find(pdgs_.begin(),pdgs_.end(),sp->pdgId());
      if(pdgfnd != pdgs_.end() ){
	auto mapfnd = emap.find(sp);
	if(mapfnd == emap.end()) 
	  emap[sp] = energy;
	else
	  emap[sp] += energy;
      }
    }
    // check if any single particle generated enough energy.  Save All the particles
    // that have enough in the map
    for(auto const& imap : emap) {
      if(debug_ > 0)std::cout << "Energy sum = " << imap.second << std::endl;
      if(imap.second >= minpe_){
	retval = true;
	output->push_back(imap.first);
      }
    }
    evt.put(std::move(output));
    // look at total energy too
    retval |= etot > minetot_;
    return retval; 
  }

}

using mu2e::CaloShowerSimFilter;
DEFINE_ART_MODULE(CaloShowerSimFilter);
