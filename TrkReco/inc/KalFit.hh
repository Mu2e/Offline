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
#include "TrkReco/inc/AmbigResolver.hh"
#include "TrkReco/inc/KalFitData.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
#include "TrkReco/inc/TrkPrintUtils.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++
#include <array>

namespace mu2e 
{
  class Calorimeter;

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
    // all functions using the KalFitData 'common block' need to be re-engineered FIXME!
// // create a fit object from a track definition
// create a fit object from  a track seed, 
    void makeTrack(StrawResponse::cptr_t srep, 
		   Mu2eDetector::cptr_t detmodel,
		   KalFitData&kalData);
// add a set of hits to an existing fit
    void addHits(StrawResponse::cptr_t srep, Mu2eDetector::cptr_t detmodel, 
		 KalFitData&kalData, double maxchi);
    // return value is the index of the cluster (if added)  
    int addTrkCaloHit(Mu2eDetector::cptr_t detmodel, KalFitData&kalData);
// add materials to a track
    bool unweedHits      (KalFitData&kalData, double maxchi);
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
    void setCalorimeter  (const Calorimeter*         Cal    ) { _calorimeter = Cal;     }
    void setTracker      (const Tracker*             Tracker) { _tracker     = Tracker; }
    void setCaloGeom();
    
    void       findCaloDiskFromTrack(KalFitData& kalData, int& trkToCaloDiskId, double&trkInCaloFlt);

    TrkErrCode fitIteration  (Mu2eDetector::cptr_t detmodel,
			      KalFitData& kalData,int iter); 
    bool       weedHits      (KalFitData& kalData, int    iter);
    bool       updateT0      (KalFitData& kalData, int    iter);
    bool       weedTrkCaloHit(KalFitData& kalData, int    iter=-1);

    bool       useTrkCaloHit() const { return _useTrkCaloHit;}
    void       fillTchDiag(KalFitData& kalData);

    bool       hit_time  (TrkHit*hit, HitT0& hitT0);
    HitT0      krep_hitT0(KalRep*krep, const TrkHit*hit);
    
    TrkPrintUtils*  printUtils() { return _printUtils; }

  private:
    // iteration-independent configuration parameters
    int _debug;		    // debug level
    double _maxhitchi;	    // maximum hit chi when adding or weeding
    double _maxpull;   // maximum pull in TrkHit 
    unsigned _maxweed;
    unsigned _maxweedtch;
    bool _initt0;	    // initialize t0?
    bool _useTrkCaloHit;    //use the TrkCaloHit 
    float  _nCaloExtrapolSteps;
    double _caloHitErr; // spatial error to use for TrkCaloHit
    std::vector<bool> _updatet0; // update t0 ieach iteration?
    std::vector<double> _t0tol;  // convergence tolerance for t0
    double _t0errfac;	    // fudge factor for the calculated t0 error
    double _mint0doca;	    // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
    double _t0nsig;	    // # of sigma to include when selecting hits for t0
    double _mindocatch, _maxdocatch; //minimum and maximum value of the TrkCaloHit DOCA
    double _mindepthtch, _maxdepthtch; //minimum and maximum value of the TrkCaloHit depth within the crystals
    double _maxtchdt; //maximum time window allowed to match a CaloCluster with the Track
    double _mintchenergy;//minimum energy of the TrkCaloHit to be considered
    double _mintchtrkpath;//minimum track path length allowed when doing the track extrapolation to the calorimeter
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
    const mu2e::Tracker*             _tracker;     // straw tracker geometry
    const mu2e::Calorimeter*         _calorimeter;
    int    _annealingStep;
    TrkTimeCalculator _ttcalc;
// relay access to BaBar field: this should come from conditions, FIXME!!!
    mutable BField* _bfield;
 
// parameters needed for evaluating the expected track impact point in the calorimeter
    unsigned _nCaloDisks;
    std::array<float,2> _zmaxcalo, _zmincalo, _rmaxcalo, _rmincalo;

    TrkPrintUtils*  _printUtils;

  // helper functions
    bool fitable(KalSeed const& kseed);
    void initT0(KalFitData&kalData);
    
    void makeTrkStrawHits  (StrawResponse::cptr_t srep, 
			    KalFitData&kalData, TrkStrawHitVector& tshv );
    void makeTrkCaloHit    (KalFitData&kalData, TrkCaloHit *&tch);
    void makeMaterials     ( Mu2eDetector::cptr_t detmodel,
			     TrkStrawHitVector const&, HelixTraj const& htraj, 
			     std::vector<DetIntersection>& dinter);
    unsigned addMaterial   (Mu2eDetector::cptr_t detmodel, KalRep* krep);
    bool unweedBestHit     (KalFitData&kalData, double maxchi);
    TrkErrCode fitTrack    (Mu2eDetector::cptr_t detmodel, KalFitData&kalData);
    void updateHitTimes    (KalRep* krep); 
    double zFlight         (KalRep* krep,double pz);
    double extendZ         (extent ex);
    TrkErrCode extendFit   (KalRep* krep);

    void findBoundingHits  (KalRep* krep, double flt0,
			    TrkHitVector::reverse_iterator& ilow,
			    TrkHitVector::iterator& ihigh);
  };
}
#endif
