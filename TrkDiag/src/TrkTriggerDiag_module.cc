//
// Stereo hit diagnostics.  Split out of MakeStereoHits
//
// Original author D. Brown
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
// root
#include "TH1F.h"
// data
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include <vector>

using namespace std; 

namespace mu2e 
{

  class TrkTriggerDiag : public art::EDAnalyzer {
    public:
      explicit TrkTriggerDiag(fhicl::ParameterSet const&);
      virtual ~TrkTriggerDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      // configuration
      int _debug, _diag;
      // diagnostics
      // helper functions
      std::vector<int> _trigbits; // list of valid trigger bits, should be in BitMap FIXME!
      // histograms
      TH1F* _mergedbits;
      std::map<std::string,TH1F*> _flagbits;
  };

  TrkTriggerDiag::TrkTriggerDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _debug		(pset.get<int>("debugLevel",0)),
    _diag		(pset.get<int>("diagLevel",1))
  {}

  TrkTriggerDiag::~TrkTriggerDiag(){}

  void TrkTriggerDiag::beginJob() {
    // create diagnostics if requested
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _mergedbits = tfs->make<TH1F>("Merged","Merged Trigger Flag;Trigger Bits",32,-0.5,31.5);
      // define which bits to test.
      _trigbits.push_back(TriggerFlag::prescaleRandom);
      _trigbits.push_back(TriggerFlag::hitCluster);
      _trigbits.push_back(TriggerFlag::helix);
      _trigbits.push_back(TriggerFlag::track);
      _trigbits.push_back(TriggerFlag::caloCluster);
    }
  }

  void TrkTriggerDiag::analyze(const art::Event& evt ) {
  // summarize all the trigger info objects
    TriggerFlag tf;
// get all the TriggerInfo objects in the event
    typedef vector< art::Handle<TriggerInfo> > HandleVector;
    HandleVector tiHV;
    evt.getManyByType(tiHV);
    for(auto tiH : tiHV) {
      if(tiH.isValid()) {
	auto const triginfo = tiH.product();
	// find the histogram associated with this module
	auto module = tiH.provenance()->moduleLabel();
	auto ifind = _flagbits.find(module);
	TH1F* bhist(0);
	if(ifind != _flagbits.end()){
	  bhist = ifind->second;
	} else {
// create a new histogram for this module
	  art::ServiceHandle<art::TFileService> tfs;
	  string title = module + string(";Trigger Bits");
	  bhist = tfs->make<TH1F>(module.c_str(),title.c_str(),32,-0.5,31.5);
	  _flagbits.insert(make_pair(module,bhist));
	}
	auto trigflag = triginfo->triggerBits();
	for(auto tbit : _trigbits) {
	  if(trigflag.hasAllProperties((TriggerFlag::bit_type)tbit))
	    bhist->Fill(tbit);
	}
	// merge this flag for the event
	tf.merge(trigflag);
	if(_debug > 0){
	  cout << "TriggerInfo for " << *tiH.provenance() 
	    <<" Has bits " << trigflag.stringRep() << endl;
	}
      }
    }
    // histogram bits of the event merger
    if(_diag > 0){
      for(auto tbit : _trigbits) {
      if(tf.hasAllProperties((TriggerFlag::bit_type)tbit))
	_mergedbits->Fill(tbit);
      }
    }
  }
}

// Part of the magic that makes this class a module.
using mu2e::TrkTriggerDiag;
DEFINE_ART_MODULE(TrkTriggerDiag);

