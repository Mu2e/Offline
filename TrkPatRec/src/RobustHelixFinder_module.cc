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
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
/// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
//#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
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
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce(art::Event& event ); 
  private:
    unsigned                           _iev;

    // configuration parameters
    int                                _diag,_debug;
    int                                _printfreq;

    art::Handle<StrawHitCollection>    _strawhitsH;
    art::Handle<TrackSeedCollection>   _trkseedsH;

    // event object labels
    std::string                        _shLabel;
    std::string                        _shpLabel;
    std::string                        _timeclusterLabel;

    // outlier cuts
    TrkParticle                        _tpart; // particle type being searched for
    TrkFitDirection                    _fdir;  // fit direction in search
    double			       _helicity; // cache the value of helicity, which can be computed from the above with the BField

    // cache of event objects
    const StrawHitCollection*          _shcol;
    const StrawHitPositionCollection*  _shpcol;
    const TimeClusterCollection*       _tccol;

    // robust helix fitter
    RobustHelixFit                     _hfit;

    // helper functions
    bool findData           (const art::Event& e);
    void fillTrackSeed      (TrackSeed &tmpseed     , 
			     TrkDef    &seeddef     ,  
			     TrackSeed  InputTrkSeed);
    
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _shLabel     (pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel    (pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _trkseedLabel(pset.get<string>("TrackSeedCollectionLabel","TimePeakFinder")),
    _tpart       ((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir        ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _helicity    (0.0),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet()))
  {
    produces<TrackSeedCollection>();
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::beginJob(){
  }

  void RobustHelixFinder::beginRun(art::Run& ){
  // calculate the helicity
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    // positive helicity is clockwise rotation around the direction of axial motion (negative dphi/dz)  
    _helicity = copysign(1.0,-_fdir.dzdt()*_tpart.charge()*field.z());
 
  }

  void RobustHelixFinder::produce(art::Event& event ) {

    // create output
    unique_ptr<TrackSeedCollection> outseeds(new TrackSeedCollection);
    // event printout
    _iev=event.id().event();
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::RobustHelixFinder: data missing or incomplete"<< endl;
    }

    for(auto tclust: _tccol) {
      // create track definitions for the helix fit from this initial information
      HelixDef       helixdef(_shcol, _shpcol, trkSeed->_timeCluster._strawHitIdxs, _tpart, _fdir );

      // copy this for the other fits
      TrkDef         seeddef(helixdef);

      // track fitting objects for this peak
      HelixFitResult helixfit(helixdef);

      // robust helix fit
      if(_hfit.findHelix(helixfit, _diag)){

	// convert the result to standard helix parameters, and initialize the seed definition helix
	HepVector hpar;
	HepVector hparerr;
	_hfit.helixParams(helixfit,hpar,hparerr);
	HepSymMatrix hcov = vT_times_v(hparerr);
	seeddef.setHelix(HelixTraj(hpar,hcov));
	// Filter outliers using this helix
	// This functionality seems to have been removed, plus it's not clear the list of hits used in the
	// helix fit is preserved when making the seed, FIXME!!
	if (_debug>1) {std::cout <<"RobustHelixFinder::produce - helix params " << hpar << "and errors " << hparerr << endl;}
	//fill seed information
	// This copying of data between related classes is error prone and unnecessary, FIXME!!
	TrackSeed tmpseed;
	fillTrackSeed(tmpseed, seeddef, *trkSeed);
	// verify the seed has the right helicity
	if(_helicity*tmpseed._helix.helicity() > 0.0){
	  outseeds->push_back(tmpseed);
	} else if (_debug > 0){
	  std::cout << "Found seed with wrong helicity " << std::endl;
	}
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

    if (evt.getByLabel(_trkseedLabel, _trkseedsH))
      _tccol = _trkseedsH.product();

    return _shcol != 0 && _shpcol != 0 && _tccol!=0;
  }


  void RobustHelixFinder::fillTrackSeed(TrackSeed &tmpseed     , 
					TrkDef    &seeddef     ,  
					TrackSeed  InputTrkSeed) {

    tmpseed._timeCluster._z0            = InputTrkSeed._timeCluster._z0;
    tmpseed._timeCluster._t0            = InputTrkSeed._timeCluster._t0;
    tmpseed._timeCluster._errt0         = InputTrkSeed._timeCluster._errt0;	  
    tmpseed._helix._d0                  = seeddef.helix().d0();
    tmpseed._helix._phi0                = seeddef.helix().phi0();
    tmpseed._helix._omega               = seeddef.helix().omega();
    tmpseed._helix._z0                  = seeddef.helix().z0();
    tmpseed._helix._tanDip              = seeddef.helix().tanDip();

    for (std::vector<hitIndex>::const_iterator ihit=seeddef.strawHitIndices().begin(); ihit!=seeddef.strawHitIndices().end(); ++ihit) {
      tmpseed._timeCluster._strawHitIdxs.push_back( mu2e::hitIndex( ihit->_index, ihit->_ambig) );
    }
  }

  void RobustHelixFinder::HelixVal2HelixTraj (const HelixVal &helIn, HelixTraj &helOut) {
    CLHEP::HepVector helParams(5);
    helParams(1) = helIn._d0;
    helParams(2) = helIn._phi0;
    helParams(3) = helIn._omega;
    helParams(4) = helIn._z0;
    helParams(5) = helIn._tanDip;
    CLHEP::HepSymMatrix conv(5,1);

    HelixTraj tmpHelix(helParams,conv);
    helOut=tmpHelix;
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
