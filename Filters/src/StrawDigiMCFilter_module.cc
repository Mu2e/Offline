//
// Select events with a minimum number of Straw Digis coming from a single
// particle, and save those particles as 'primaries'.  Particle can also be
// selected based on PDG code and momentum
//  original author: David Brown (LBNL)
//

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class StrawDigiMCFilter : public art::EDFilter {
    public:
      explicit StrawDigiMCFilter(fhicl::ParameterSet const& pset);
      virtual ~StrawDigiMCFilter() { }

      bool filter( art::Event& event);

    private:

      unsigned minndigi_;
      double minpmom_, maxpmom_;
      size_t minnparticles_;
      std::vector<PDGCode::type> pdgs_;
      int diag_, debug_;
      art::InputTag _mcdigisTag;

  };

  StrawDigiMCFilter::StrawDigiMCFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    minndigi_(pset.get<unsigned>("MinNDigis")),
    minpmom_(pset.get<double>("MinParticleMom")),
    maxpmom_(pset.get<double>("MaxParticleMom")),
    minnparticles_(pset.get<size_t>("MinNParticles", 1)),
    diag_(pset.get<int>("diagLevel",0)),
    debug_(pset.get<int>("debugLevel",0)),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")){
      produces<SimParticlePtrCollection>();
      auto pdgs = pset.get<std::vector<int>>("particleTypes");
      for(auto pdg : pdgs )
        pdgs_.push_back(PDGCode::type(pdg));
    }


  bool StrawDigiMCFilter::filter(art::Event& evt) {
    std::unique_ptr<SimParticlePtrCollection> output(new SimParticlePtrCollection());
    bool retval(false);
    auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
    const StrawDigiMCCollection *mcdigis = mcdH.product();
    // keep count of digis produced by specific particle
    std::map<art::Ptr<SimParticle>,unsigned> pmap;
    for(auto const& mcdigi : *mcdigis) {
      // do not inspect truth info for digis not produced in simulation
      if (mcdigi.provenance().ContainsSimulation()){
        // look at the early end
        StrawEnd fend = mcdigi.earlyEnd();
        auto const& step =  mcdigi.strawGasStep(fend);
        art::Ptr<SimParticle> const& sp = step->simParticle();
        auto const& mom = step->momentum(); // cast to 3-vector
        if(debug_ > 0)std::cout <<"SimParticle PDG = " << sp->pdgId() << " Mom = " << sqrt(mom.mag2()) << std::endl;
        bool goodpdg(true);
        if(pdgs_.size() > 0)
          goodpdg = std::find(pdgs_.begin(),pdgs_.end(),sp->pdgId()) != pdgs_.end();
        if(goodpdg && sqrt(mom.mag2()) > minpmom_ && sqrt(mom.mag2()) < maxpmom_ ){
          auto mapfnd = pmap.find(sp);
          if(mapfnd == pmap.end())
            pmap[sp] = 1;
          else
            ++mapfnd->second;
        }
      }
    }
    // check if any single particle generated enough digis.  Save All the particles
    // that have enough in the map
    size_t nParticles(0);
    for(auto const& imap : pmap) {
      if(imap.second >= minndigi_){
        retval = true;
        output->push_back(imap.first);
        ++nParticles;
      }
    }
    if (nParticles < minnparticles_) retval = false;
    evt.put(std::move(output));
    return retval;
  }

}

using mu2e::StrawDigiMCFilter;
DEFINE_ART_MODULE(StrawDigiMCFilter)
