//
// Class to perform BaBar Kalman fit
// Original author: Dave Brown LBNL 2012
//
// $Id: KalFit.hh,v 1.29 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef TrkReco_KalFit_HH
#define TrkReco_KalFit_HH

// framework

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
// Mu2e objects
#include "Mu2eBTrk/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "TrkReco/inc/KalFitResult.hh"
#include "TrkReco/inc/AmbigResolver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++

namespace mu2e 
{
  class KalFit : public KalContext
  {
  public:
// define different ambiguity resolution strategies
    enum ambigStrategy {fixedambig=0,pocaambig=1,hitambig=2,panelambig=3,doubletambig=4};
// parameter set should be passed in on construction
#ifndef __GCCXML__
    explicit KalFit(fhicl::ParameterSet const&, TrkDirection const& tdir);
#endif/*__GCCXML__*/

    virtual ~KalFit();
// main function: given a track definition, create a fit object from it
    virtual void makeTrack(KalFitResult& kres);
// add a set of hits to an existing fit
    virtual void addHits(KalFitResult& kres,const StrawHitCollection* straws, std::vector<hitIndex> indices, double maxchi);
// Try to put back inactive hits
    bool unweedHits(KalFitResult& kres, double maxchi);
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
  private:
    // iteration-independent configuration parameters
    int _debug;		    // debug level
    double _maxhitchi;	    // maximum hit chi when adding or weeding
    unsigned _maxweed;	    // maximum number of hits to remove when weeding
    double _maxdriftpull;   // maximum drift pull in TrkStrawHit 
    bool _initt0;	    // initialize t0?
    bool _updatet0;	    // update t0 ieach iteration?
    double _t0errfac;	    // fudge factor for the calculated t0 error
    double _mint0doca;	    // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
    double _t0nsig;	    // # of sigma to include when selecting hits for t0
    bool _removefailed;	    // remove failed fits?
    unsigned _minnstraws;   // minimum # staws for fit
    // iteration-dependent configuration parameters
    std::vector<bool> _weedhits;	// weed hits?
    std::vector<double> _herr;		// what external hit error to add (for simulated annealing)
    std::vector<int> _ambigstrategy;	// which ambiguity resolver to use
    std::vector<AmbigResolver*> _ambigresolver;	
    std::vector<double> _t0tol;  // convergence tolerance for t0
// state
    TrkParticle _tpart;
    TrkFitDirection _fdir;
// relay access to BaBar field: this should come from conditions, FIXME!!!
    mutable BField* _bfield;
  // helper functions
    bool fitable(TrkDef const& tdef);
    void initT0(TrkDef const& tdef, TrkT0& t0);
    
    bool weedHits(KalFitResult& kres); //  KalRep and TrkHit
    bool updateT0(KalFitResult& kres); // KalRep and TrkStrawHit
    void fitTrack(KalFitResult& kres); // KalRep
    virtual void makeHits(KalFitResult& kres, TrkT0 const& t0); // TrkDef and TrkStrawHit
    virtual void makeMaterials(KalFitResult& kres); // TrkDef and TrkStrawHit
    void fitIteration(KalFitResult& kres); // KalRep and TrkStrawHit
    void updateHitTimes(KalFitResult& kres); // KalRep and TrkStrawHit
    
    void findBoundingHits(std::vector<TrkStrawHit*>& hits, double flt0,
	std::vector<TrkStrawHit*>::reverse_iterator& ilow,
	std::vector<TrkStrawHit*>::iterator& ihigh);
  };
}
#endif
