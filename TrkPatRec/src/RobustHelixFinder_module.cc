//
// TTracker Pattern Recognition based on Robust Helix Fit
//
// $Id: RobustHelixFinder_module.cc,v 1.2 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
// conditions
/// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
using namespace std; 
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e 
{
  class RobustHelixFinder : public art::EDProducer
  {
  public:
    explicit RobustHelixFinder(fhicl::ParameterSet const&);
    virtual ~RobustHelixFinder();
    virtual void produce(art::Event& event ); 
  private:
    unsigned                           _iev;

    // configuration parameters
    int                                _diag,_debug;
    int                                _printfreq;
    bool				_saveall;

    // input object tags
    art::InputTag			_shTag;
    art::InputTag			_shpTag;
    art::InputTag			_shfTag;
    art::InputTag			_tcTag;
    // output label
    std::string                        _trackseed;

    // hit selection
    StrawHitFlag  _psel;

    // cache of event objects
    const StrawHitCollection*          _shcol;
    const StrawHitPositionCollection*  _shpcol;
    const StrawHitFlagCollection*      _shfcol;
    const TimeClusterCollection*       _tccol;

    // robust helix fitter
    RobustHelixFit                     _hfit;

    // helper functions
    bool findData           (const art::Event& e);
    
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _saveall     (pset.get<bool>("SaveAllHelices",false)),
    _shTag	 (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	 (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	 (pset.get<art::InputTag>("StrawHitFlagCollection","TimeClusterFinder")),
    _tcTag	 (pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _trackseed   (pset.get<string>("HelixSeedCollectionLabel","TimeClusterFinder")),
    _psel        (pset.get<std::vector<std::string> >("PositionSelectionBits")),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet()))
  {
    produces<HelixSeedCollection>();
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::produce(art::Event& event ) {

    // create output
    unique_ptr<HelixSeedCollection> outseeds(new HelixSeedCollection);
    // event printout
    _iev=event.id().event();
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::RobustHelixFinder: data missing or incomplete"<< endl;
    }

    for(auto const& tclust: *_tccol) {
    // build an empty HelixSeed 
    // loop over hits in this time cluster and select  hits with good 3-d position information
    std::vector<hitIndex> goodhits;
      for(auto const& ind : tclust._strawHitIdxs) {
	if(_shfcol->at(ind).hasAnyProperty(_psel))
	  goodhits.push_back(ind);
      }
     // build a helix seed using these hits, but the original t0
      HelixSeed hseed;
      hseed._timeCluster = tclust;
      hseed._timeCluster._strawHitIdxs = goodhits;
      _hfit.findHelix(*_shcol, *_shpcol, hseed);

// should iterate fit to include outlier removal using time + geometric information FIXME!

      if((hseed._status.hasAllProperties(TrkFitFlag::fitOK) && _hfit.helicity() == hseed._helix.helicity()) || _saveall) {
	  outseeds->push_back(hseed);
	  if(_debug > 1) cout << "Found helix with fit \n" << hseed._helix << endl;
      } else if (_debug > 1) cout << "Found helix without fit \n" << hseed._helix << endl;
    }

    if (_debug>0 && (_iev%_printfreq)==0) {
      cout<<"event "<<_iev<<" tot N hit "<<_shcol->size()<<" N tracks seed found "<<outseeds->size()
	       <<" N time peaks "<<_tccol->size()<<endl;
    }

    event.put(std::move(outseeds));
  }

  // find the input data objects 
  bool RobustHelixFinder::findData(const art::Event& evt){
    _shcol = 0; _shfcol = 0; _shpcol = 0; _tccol = 0; 
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    _tccol = tcH.product();

    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && _tccol != 0;
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
