// Select events with
// 1) a minimum number of Straw Digis coming from a single particle if it has a momentum > minimum momentum
// 2) a minimum number of "hit" stations if the momentum of the particle < minimum momentum
//  Adapted from the CosmicStrawDigiMCFilter module
//

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
// data
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class CosmicHelixMCFilter : public art::EDFilter {
    public:
      explicit CosmicHelixMCFilter(fhicl::ParameterSet const& pset);
      virtual ~CosmicHelixMCFilter() { }

      bool filter(art::Event& event);

    private:

    int _diag, _debug,_checkpid;
    //art::ProductToken<ComboHitCollection> const _shToken;
    art::InputTag _hsToken;
    // cache of event objects
    //const ComboHitCollection *_chcol;
    const HelixSeedCollection *_hscol; 
    art::InputTag _mcDigisTag;
  
  };

  CosmicHelixMCFilter::CosmicHelixMCFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _diag(pset.get<int>("diag",0)),
    _debug(pset.get<int>("debug",0)),
    _checkpid(pset.get<int>("checkpid",13)),
    //_shToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    //_hsToken{consumes<HelixSeedCollection>(pset.get<art::InputTag>("SeedCollection"))},
    /*{ consumes<HelixSeedCollection>(_hsToken);
      consumes<StrawDigiMCCollection>(_mcDigisTag);
      }*/
    _hsToken(pset.get<art::InputTag>("HelixSeedCollection","MHFinderDmu")),
    _mcDigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")){
      produces<HelixSeedCollection>();
    }

    //_mcDigisTag{consumes<StrawDigiMCCollection>(pset.get<art::InputTag>("StrawDigiMCCollection"))}
    
  bool CosmicHelixMCFilter::filter(art::Event& evt) {
    bool retval(false);
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsToken);
    _hscol = hsH.product();
    if(_debug > 0)std::cout<<" helix size = "<<_hscol->size()<<std::endl;
    auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcDigisTag);
    const StrawDigiMCCollection *mcdigis = mcdH.product();
    std::unique_ptr<HelixSeedCollection> output(new HelixSeedCollection());
    // loop over the Helices
    if(_hscol->size()>=1){
      for (size_t iseed=0; iseed<_hscol->size(); ++iseed) {
        // convert the HelixSeed to a TrkDef
        HelixSeed const& hseed(_hscol->at(iseed));
        if(_debug>0)std::cout<<"Helix combo hit size = "<<hseed.hits().size()<<std::endl;
        //int nss_ch = hseed.hits().size();
        int loc, pdgid(0);
        for(uint16_t ihit=0;ihit < hseed.hits().size(); ++ihit){
          ComboHit const& ch = hseed.hits()[ihit];
          loc   = ch.index();
	  //if(_debug>0) std::cout<<"LOC = "<<loc<<std::endl;
          if (loc >= 0) {
	      auto const& mcdigi = mcdigis->at(loc);
              StrawEnd fend = mcdigi.earlyEnd();
              auto const& step =  mcdigi.strawGasStep(fend);
              art::Ptr<SimParticle> const& sp = step->simParticle();
              int pid = sp->pdgId();
	      pdgid = pdgid + fabs(pid);
	  }
	}
	if(pdgid >0) pdgid = pdgid/hseed.hits().size();
	if(_debug >0) std::cout<<"Final PID of helix = "<<pdgid<<std::endl;
        if(pdgid == _checkpid)  {output->push_back(hseed);
            retval = true;
	  }
      }
    }
    if(output->size()>0 and _debug>0) std::cout<<"output size = "<<output->size()<<std::endl;
    evt.put(std::move(output));
    return retval;
  }

}

using mu2e::CosmicHelixMCFilter;
DEFINE_ART_MODULE(CosmicHelixMCFilter)
