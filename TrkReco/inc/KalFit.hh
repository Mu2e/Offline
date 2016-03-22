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
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
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
    enum ambigStrategy {fixedambig=0,pocaambig=1,hitambig=2,panelambig=3,doubletambig=4};
// different locations to which the track may be extended
    enum extent {noextension=-1,target=0,ipa=1,tracker=2,calo=3};
// parameter set should be passed in on construction
#ifndef __GCCXML__
    explicit KalFit(fhicl::ParameterSet const&, TrkFitDirection const& tdir);
#endif/*__GCCXML__*/

    virtual ~KalFit();
// main function: given a track definition, create a fit object from it
    virtual void makeTrack(TrkDef const& tdef, KalRep*& kres);
// add a set of hits to an existing fit
    virtual void addHits(KalRep* kres,const StrawHitCollection* straws, std::vector<hitIndex> indices, double maxchi);
// add materials to a track
    bool unweedHits(KalRep* kres, double maxchi);
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
  private:
    // iteration-independent configuration parameters
    int _debug;		    // debug level
    double _maxhitchi;	    // maximum hit chi when adding or weeding
    double _maxdriftpull;   // maximum drift pull in TrkStrawHit 
    bool _initt0;	    // initialize t0?
    bool _updatet0;	    // update t0 ieach iteration?
    std::vector<double> _t0tol;  // convergence tolerance for t0
    double _t0errfac;	    // fudge factor for the calculated t0 error
    double _mint0doca;	    // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
    double _t0nsig;	    // # of sigma to include when selecting hits for t0
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
// state
    TrkParticle _tpart;
    TrkFitDirection _fdir;
// relay access to BaBar field: this should come from conditions, FIXME!!!
    mutable BField* _bfield;
  // helper functions
    bool fitable(TrkDef const& tdef);
    void initT0(TrkDef const& tdef, TrkT0& t0);
    virtual void makeHits(TrkDef const& tdef, TrkT0 const& t0, TrkStrawHitVector& tshv); 
    virtual void makeMaterials(TrkStrawHitVector const&, TrkDef const& tdef, std::vector<DetIntersection>& dinter);
    unsigned addMaterial(KalRep* krep);
    bool weedHits(KalRep* kres, TrkStrawHitVector& tshv,size_t iter);
    bool unweedHits(KalRep* kres, TrkStrawHitVector& tshv, double maxchi);
    bool updateT0(KalRep* kres, TrkStrawHitVector& tshv);
    TrkErrCode fitTrack(KalRep* kres, TrkStrawHitVector& tshv);
    TrkErrCode fitIteration(KalRep* kres,TrkStrawHitVector& tshv,size_t iter); 
    void updateHitTimes(KalRep* kres, TrkStrawHitVector& tshv); 
    double zFlight(KalRep* krep,double pz);
    double extendZ(extent ex);
    TrkErrCode extendFit(KalRep* krep);

    void findBoundingHits(TrkStrawHitVector& hits, double flt0,
	TrkStrawHitVector::reverse_iterator& ilow,
	TrkStrawHitVector::iterator& ihigh);
  };
}
#endif
