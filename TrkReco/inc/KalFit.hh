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
#include "RecoDataProducts/inc/ComboHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
// Mu2e objects
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrkData/inc/TrkCaloHit.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "TrkReco/inc/TrkDef.hh"
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
    enum ambigStrategy {fixedambig=0,hitambig=2,panelambig=3,doubletambig=4};
// different locations to which the track may be extended
    enum extent {noextension=-1,target=0,ipa=1,tracker=2,calo=3};
// parameter set should be passed in on construction
#ifndef __GCCXML__
    explicit KalFit(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/

    virtual ~KalFit();
// // create a fit object from a track definition
//     void makeTrack(const ComboHitCollection* shcol, TrkDef& tdef, KalRep*& kres);
// create a fit object from  a track seed, 
    void makeTrack(const ComboHitCollection* shcol, KalSeed const& kseed, KalRep*& kres);
// add a set of hits to an existing fit
    void addHits(KalRep* kres,const ComboHitCollection* shcol, std::vector<StrawHitIndex> indices, double maxchi);
// add materials to a track
    bool unweedHits(KalRep* kres, double maxchi);
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
  private:
    // iteration-independent configuration parameters
    int _debug;		    // debug level
    double _maxhitchi;	    // maximum hit chi when adding or weeding
    double _maxpull;   // maximum pull in TrkHit 
    bool _initt0;	    // initialize t0?
    bool _useTrkCaloHit;    //use the TrkCaloHit to initialize the t0?
    bool _updatet0;	    // update t0 ieach iteration?
    std::vector<double> _t0tol;  // convergence tolerance for t0
    double _t0errfac;	    // fudge factor for the calculated t0 error
    double _mint0doca;	    // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
    double _t0nsig;	    // # of sigma to include when selecting hits for t0
    double _dtoffset;
    double _strHitW, _calHitW;//weight used to evaluate the initial track T0
    unsigned _minnstraws;   // minimum # staws for fit
    double _maxmatfltdiff; // maximum difference in track flightlength to separate to intersections of the same material
    // iteration-dependent configuration parameters
    std::vector<bool> _weedhits;	// weed hits?
    std::vector<double> _herr;		// what external hit error to add (for simulated annealing)
    std::vector<int> _ambigstrategy;	// which ambiguity resolver to use
    std::vector<bool> _addmaterial; // look for additional materials along the track
    std::vector<AmbigResolver*> _ambigresolver;
    bool _resolveAfterWeeding;
    extent _exup;
    extent _exdown;
// relay access to BaBar field: this should come from conditions, FIXME!!!
    mutable BField* _bfield;
  // helper functions
    bool fitable(TrkDef const& tdef);
    bool fitable(KalSeed const& kseed);
    void initT0(KalRep* krep);
    
    void makeTrkStrawHits  (const ComboHitCollection* shcol, HelixTraj const& htraj,
			    std::vector<TrkStrawHitSeed>const& hseeds, TrkStrawHitVector& tshv );
    void makeTrkCaloHit    (KalSeed const& kseed, TrkCaloHit *&tch);
    void makeMaterials     (TrkStrawHitVector const&, HelixTraj const& htraj, std::vector<DetIntersection>& dinter);
    unsigned addMaterial   (KalRep* krep);
    bool weedHits          (KalRep* kres, size_t iter);
    bool unweedBestHit     (KalRep* kres, double maxchi);
    bool updateT0          (KalRep* kres);
    TrkErrCode fitTrack    (KalRep* kres);
    TrkErrCode fitIteration(KalRep* kres,size_t iter); 
    void updateHitTimes    (KalRep* kres); 
    double zFlight         (KalRep* krep,double pz);
    double extendZ         (extent ex);
    TrkErrCode extendFit   (KalRep* krep);

    void findBoundingHits  (KalRep* krep, double flt0,
			    TrkHitVector::reverse_iterator& ilow,
			    TrkHitVector::iterator& ihigh);
    
    void findTrkCaloHit    (KalRep*krep, TrkCaloHit*tch);
  };
}
#endif
