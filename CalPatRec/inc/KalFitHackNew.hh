//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitHackNew.hh,v 1.3 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
//
#ifndef KalFitHackNew_HH
#define KalFitHackNew_HH

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
// tracker
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
// KalFitHackNew objects
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

#include "TrkReco/inc/AmbigResolver.hh"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/Doublet.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"

//ROOT

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {

  class Calorimeter;
  class KalFitHackNew : public KalContext {
  public:
					// different ambiguity resolution strategies
    enum ambigStrategy {
      kFixedAmbig    = 0,
      kPocaAmbig     = 1,
      kHitAmbig      = 2,
      kPanelAmbig    = 3,
      kDoubletAmbig  = 4
    };
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
  protected:
    // configuration parameters
    int             

            _debugLevel;
    unsigned                    _minnstraws;
    double                      _maxmatfltdiff; // maximum difference in track flightlength to separate to intersections of the same material
    vector<bool>                _weedhits;
    double                      _maxhitchi;
    unsigned                    _maxweed;
    std::vector<double>         _hiterr;
    double                      _maxdriftpull;
    fhicl::ParameterSet*        _darPset;         // parameter set for doublet ambig resolver
    std::vector<AmbigResolver*> _ambigresolver;
    //    bool                        _initt0;
    bool                        _updateT0;
    int                         _updateT0Mode;    // 0: use cluster T0 1:update T0 assuming no cluster time
    double                      _minHitDrift;
    double                      fRdriftMinusDocaTol;
    int                         fSign[4][2];
    std::vector<double>         _t0tol;
    double                      _t0errfac;       // fudge factor for the calculated t0 error
    double                      _mint0doca;      // minimum (?) doca for t0 hits
    double                      _t0nsig;         // # of sigma to include when selecting hits for t0
    double                      _dtoffset;       // track - luster time offset, ns
    double                      fScaleErrDoublet;
    double                      fMinDriftDoublet;
    double                      fDeltaDriftDoublet;
    double                      _maxDoubletChi2;
    double                      fSigmaSlope;
    std::string                 fMakeStrawHitModuleLabel;
		                
    bool                        _removefailed;
    TrkParticle                 _tpart;
    TrkFitDirection             _fdir;
    std::vector<int>            _ambigstrategy;
    std::vector<bool>           _addmaterial;
    mutable BField*             _bfield;
    int                         _nIter;
    const CalTimePeak*          fTimePeak;
    int                         _annealingStep;

    const mu2e::PtrStepPointMCVectorCollection*  _listOfMCStrawHits;

    const mu2e::Tracker*             _tracker;     // straw tracker geometry
    const mu2e::TrackerCalibrations* _tcal;

    const mu2e::Calorimeter*         _calorimeter;
//-----------------------------------------------------------------------------
// to decouple from MC classes, need a redefinable function here
//-----------------------------------------------------------------------------
    double (*_MCDoca) (const art::Event* Event, const char* ShDigiLabel, const Straw* aStraw);
//-----------------------------------------------------------------------------
// constructors and destructor, parameter set should be passed in on construction
//-----------------------------------------------------------------------------
  public:
    explicit KalFitHackNew(fhicl::ParameterSet const&);
    virtual ~KalFitHackNew();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
    int  nIter       () { return _nIter        ; } 
    int  maxIteration() { return _hiterr.size()-1; }
    int  minNStraws  () { return _minnstraws;      }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
    void setNIter        (int N   ) { _nIter        = N  ; }

    void setTracker      (const Tracker*             Tracker) { _tracker     = Tracker; }
    void setTrackerCalib (const TrackerCalibrations* TCal   ) { _tcal        = TCal;    }

    void setCalorimeter  (const Calorimeter*         Cal    ) { _calorimeter = Cal;     }

    void setStepPointMCVectorCollection(const mu2e::PtrStepPointMCVectorCollection* List) {
      _listOfMCStrawHits = List;
    }
//-----------------------------------------------------------------------------
// add a set of hits to an existing fit
//-----------------------------------------------------------------------------
    virtual void addHits(KalFitResultNew& KRes, double MaxChi);

    void findBoundingHits(std::vector<TrkHit*>&                   hits, 
			  double                                  flt0,
			  std::vector<TrkHit*>::reverse_iterator& ilow ,
			  std::vector<TrkHit*>::iterator&         ihigh);

    //    bool fitable     (KalFitResultNew& Kres);
    void fitIteration(KalFitResultNew& Kres, int Iteration);
    void fitTrack    (KalFitResultNew& Kres);
    //    void initCaloT0  (KalFitResultNew& KRes, const CaloCluster* Cluster);
    void initT0      (KalFitResultNew& KRes);
//-----------------------------------------------------------------------------
// main function: given a track definition, create a fit object from it
//-----------------------------------------------------------------------------
    virtual void makeTrack      (KalFitResultNew& kRes);
//---------------------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following function for printing the hits included in the track
//----------------------------------------------------------------------------------------------
    void printHits (KalFitResultNew& kres, const char* Message);  // default

    bool unweedHits       (KalFitResultNew& kres, double maxchi);

    void updateCalT0      (KalFitResultNew& KRes);
    bool updateT0         (KalFitResultNew& KRes);

    bool weedHits         (KalFitResultNew& KRes, int Iteration);

// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
//-----------------------------------------------------------------------------
// functions taking KalRep* , not KalFitResultNew
//-----------------------------------------------------------------------------
    void updateHitTimes  (KalRep* KRep);
    void makeHits        (KalFitResultNew& KRes, TrkHitVector& ListOfHits);

    void makeMaterials   (TrkHitVector&            ListOfHits, 
			  const HelixTraj&         Helix, 
			  vector<DetIntersection>& ListOfIntersections);
//-----------------------------------------------------------------------------
// latest code by Dave
//-----------------------------------------------------------------------------
    double   zFlight    (KalRep* KRep, double pz); // *** FIXME category
    unsigned addMaterial(KalRep* KRep);

//-----------------------------------------------------------------------------
// default plugin
//-----------------------------------------------------------------------------
    static double MCDoca(const art::Event* Event, const char* ShDigiLabel, const Straw* aStraw) { return -99.; }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
    void setMCDocaRoutine(void* Function) { 
      _MCDoca = (double (*) (const art::Event* Event, const char* ShDigiLabel, const Straw* aStraw)) Function; 
    }
  };
}
#endif
