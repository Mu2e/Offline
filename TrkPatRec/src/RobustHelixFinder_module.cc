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
using CLHEP::Hep3Vector;
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
    unsigned				_maxniter;
    double				_cradres; // average center resolution along center position (mm)
    double				_cperpres; // average center resolution perp to center position (mm)
    double				_radres; // average radial resolution for circle fit (mm)
    double				_phires; // average azimuthal resolution on circle (rad)
    double				_maxdwire; // outlier cut on distance between hit and helix along wire
    double				_maxdtrans; // outlier cut on distance between hit and helix perp to wire
    double				_maxchisq; // outer cut on chisquared

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
    const StrawHitCollection*  _shcol;
    const StrawHitPositionCollection*  _shpcol;
    const StrawHitFlagCollection*      _shfcol;
    const TimeClusterCollection*       _tccol;

    // robust helix fitter
    RobustHelixFit                     _hfit;

    // helper functions
    bool findData           (const art::Event& e);
    bool filterHits(HelixSeed& hseed); // return value tells if any hits changed state
    void updateT0(HelixSeed& hseed); // update T0 value based on current good hits
 
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _saveall     (pset.get<bool>("SaveAllHelices",false)),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
    _cradres	 (pset.get<double>("CenterRadialResolution",10.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",10.0)),
    _radres	 (pset.get<double>("RadiusResolution",10.0)),
    _phires	 (pset.get<double>("AzimuthREsolution",0.1)),
    _maxdwire    (pset.get<double>("MaxWireDistance",100.0)), // max distance along wire
    _maxdtrans   (pset.get<double>("MaxTransDistance",100.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<double>("MaxChisquared",100.0)), // max chisquared
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
      HelixSeed hseed;
      // copy in the t0 and cluster
      hseed._t0 = tclust._t0;
      hseed._caloCluster = tclust._caloCluster;
    // loop over hits in this time cluster and select  hits with good 3-d position information
      std::vector<StrawHitIndex> goodhits;
      for(auto const& ind : tclust._strawHitIdxs) {
	if(_shfcol->at(ind).hasAnyProperty(_psel))
	  goodhits.push_back(ind);
      }
      // create helix seed hits from the straw hit positions 
      for(auto idx : goodhits ) {
	HelixHit hhit(_shpcol->at(idx),idx);
	hhit._flag.clear(StrawHitFlag::resolvedphi);
	hseed._hhits.push_back(hhit);
      }
      // iteratively fit the helix and filter outliers in space and time
      unsigned niter(0);
      bool changed(true);
      do {
	++niter;
	_hfit.findHelix(hseed);
	changed =filterHits(hseed);
	if(changed)updateT0(hseed);
      } while(hseed._status.hasAllProperties(TrkFitFlag::helixOK)  && niter < _maxniter && changed);
      // final test
      if((hseed._status.hasAllProperties(TrkFitFlag::helixOK) && _hfit.helicity() == hseed._helix.helicity()) || _saveall) {
	outseeds->push_back(hseed);
	if(_debug > 1) cout << "Found helix with fit \n" << hseed._helix << endl;
      } else if (_debug > 1) cout << "Found helix without fit \n" << hseed._helix << endl;
    }

    if (_debug>0 && (_iev%_printfreq)==0) {
      cout<<"event "<<_iev<<" tot N hit "<<_shfcol->size()<<" N tracks seed found "<<outseeds->size()
	       <<" N time peaks "<<_tccol->size()<<endl;
    }

    event.put(std::move(outseeds));
  }

  bool RobustHelixFinder::filterHits(HelixSeed& hseed) {
    RobustHelix& helix = hseed. _helix;
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
// loop over hits
    for(auto& hhit : hhits) {
      bool oldout = !hhit._flag.hasAnyProperty(outlier);
      // compute spatial distance and chisquiared
      // first, find the expected helix position for this hit's z position
      Hep3Vector hpos;
      hpos.setZ(hhit.pos().z());
      helix.position(hpos);
      Hep3Vector dh = hhit.pos() - hpos;
      double dwire = dh.dot(hhit.wdir()); // projection along wire direction
      Hep3Vector wtdir = zaxis.cross(hhit.wdir()); // transverse direction to the wire
      Hep3Vector cdir = (hhit.pos() - helix.center()).perpPart().unit(); // direction from the circle center to the hit
      Hep3Vector cperp = zaxis.cross(cdir); // direction perp to the radius
      double dtrans = dh.dot(wtdir); // transverse projection
      // compute the total resolution including hit and helix parameters
      double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	std::pow(_cradres*cdir.dot(hhit.wdir()),(int)2) +
	std::pow(_cperpres*cperp.dot(hhit.wdir()),(int)2);
      double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wtdir),(int)2);
	// need to add uncertainty in azimuth FIXME!!
      double chisq = sqrt( dwire*dwire/wres2 + dtrans*dtrans/wtres2 );
//      double dt = _shcol->at(hhit._shidx).time() - hseed._t0.t0();
      // need to add time information FIXME!
      // need to put all these variables in an MVA FIXME!!
      if(fabs(dwire) > _maxdwire || fabs(dtrans) > _maxdtrans || chisq > _maxchisq) {
	// outlier hit flag it
	hhit._flag.merge(outlier);
      } else {
	// clear the hit in case it was formerly an outlier
	hhit._flag.clear(outlier);
      }
      changed |= oldout != hhit._flag.hasAnyProperty(outlier);
 
    }

    return changed;
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

  void RobustHelixFinder::updateT0(HelixSeed& hseed) {
    // need an explicit t0 finder here to update T0 FIXME!
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
