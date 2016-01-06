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
//#include "art/Framework/Services/Optional/TFileService.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
//#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "TrkReco/inc/TrkStrawHit.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
// Mu2e
#include "TrkPatRec/inc/TrkPatRecUtils.hh"
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
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
    private:
      unsigned _iev;
      // configuration parameters
      int _diag,_debug;
      int _printfreq;
      art::Handle<mu2e::StrawHitCollection> _strawhitsH;
      art::Handle<TrackerHitTimeClusterCollection> _tclusthitH;
      // event object labels
      std::string _shLabel;
      std::string _shpLabel;
      std::string _tpkfLabel;
      // outlier cuts
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const TrackerHitTimeClusterCollection* _tccol;
      // robust helix fitter
      RobustHelixFit _hfit;
      // cache of time peaks
      std::vector<TrkTimePeak> _tpeaks;
      //
      // helper functions
      bool findData(const art::Event& e);
      Int_t _eventid;
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _tpkfLabel(pset.get<string>("TrackerHitTimeClusterCollection","TimePeakFinder")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _hfit(pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet()))
  {
     produces<TrackSeedCollection>();
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::beginJob(){
    // create a histogram of throughput: this is a basic diagnostic that should ALWAYS be on
    //art::ServiceHandle<art::TFileService> tfs;
    _eventid = 0;
  }

  void RobustHelixFinder::beginRun(art::Run& ){}

  void RobustHelixFinder::produce(art::Event& event ) {
    _eventid = event.event();

    // create output
    unique_ptr<TrackSeedCollection> outseeds(new TrackSeedCollection);
    // event printout
    _iev=event.id().event();
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::RobustHelixFinder: data missing or incomplete"<< endl;
    }
    loadTimePeaks(_tpeaks,_tccol);

    // dummy objects
    static HelixDef dummyhdef;
    static HelixFitResult dummyhfit(dummyhdef);
    for(unsigned ipeak=0;ipeak<_tpeaks.size();++ipeak){
      // create track definitions for the helix fit from this initial information 
      HelixDef helixdef(_shcol,_shpcol,_tpeaks[ipeak]._trkptrs,_tpart,_fdir);
      // copy this for the other fits
      TrkDef seeddef(helixdef);
      // track fitting objects for this peak
      HelixFitResult helixfit(helixdef);
      // robust helix fit
      if(_hfit.findHelix(helixfit/*,_icepeak==(int)ipeak*/)){
	//findhelix = true;
	// convert the result to standard helix parameters, and initialize the seed definition helix
	HepVector hpar;
	HepVector hparerr;
	_hfit.helixParams(helixfit,hpar,hparerr);
	HepSymMatrix hcov = vT_times_v(hparerr);
	seeddef.setHelix(HelixTraj(hpar,hcov));
	// Filter outliers using this helix
	if (_debug>1) {std::cout <<"RobustHelixFinder::produce - helix params " << hpar << "and errors " << hparerr << endl;}
	//fill seed information
	TrackSeed tmpseed;
	fillTrackSeed(tmpseed, seeddef, _tclusthitH, ipeak, _strawhitsH);
        outseeds->push_back(tmpseed);
      }
    }

    if (_debug>0 && (_iev%_printfreq)==0) {
            std::cout<<"event "<<_iev<<" tot N hit "<<_shcol->size()<<" N tracks seed found "<<outseeds->size()
                            <<" N time peaks "<<_tccol->size()<<std::endl;
    }

    event.put(std::move(outseeds));
  }

  // find the input data objects 
  bool RobustHelixFinder::findData(const art::Event& evt){
    _shcol = 0;
    _shpcol = 0;
    _tccol = 0;

    if(evt.getByLabel(_shLabel,_strawhitsH))
      _shcol = _strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
      _shpcol = shposH.product();

    if (evt.getByLabel(_tpkfLabel,_tclusthitH))
      _tccol = _tclusthitH.product();
// don't require stereo hits: they are only used for diagnostics
    return _shcol != 0 && _shpcol != 0 && _tccol!=0;
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
