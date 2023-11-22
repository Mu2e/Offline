// Purpose: Event filter for pair converions that make it through the tracker
// author: S Middleton and H Jafree, 2023
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class PairConvFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> StrawToken{Name("StrawGasStepCollection"),Comment("tag for straw gas step collection")};
      };
      explicit PairConvFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;
      //virtual bool beginRun(art::Run&   run   );
      //virtual bool endRun( art::Run& run ) override;
    private:
	  	art::InputTag _StrawToken;
	  	const StrawGasStepCollection* _StrawCol;
  };

  PairConvFilter::PairConvFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config},
   _StrawToken(config().StrawToken())
  {}

  bool PairConvFilter::filter(art::Event& event) {
    
    bool pass = false;
          //------------ Straw Gas Collection  --------------------//
	auto strawH = event.getValidHandle<StrawGasStepCollection>(_StrawToken);
	_StrawCol = strawH.product();
	unsigned int nElectronSteps = 0;
	unsigned int nPositronSteps = 0;
	for (size_t k = 0; k < _StrawCol->size(); k++){
  	StrawGasStep strawgas = (*_StrawCol)[k];
  	art::Ptr<SimParticle> const& simpart = strawgas.simParticle();
  		if      (simpart->pdgId()==11)  ++nElectronSteps;
  		else if (simpart->pdgId()==-11) ++nPositronSteps;
		}      

    if (nElectronSteps>=10 and nPositronSteps>=10) pass=true;
    return pass;
  }
}

using mu2e::PairConvFilter;
DEFINE_ART_MODULE(PairConvFilter)
