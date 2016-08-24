//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitHack.hh,v 1.3 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
//
#ifndef KalFitHack_HH
#define KalFitHack_HH

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
// tracker
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
// KalFitHack objects
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
#include "CalPatRec/inc/KalFitResult.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "RecoDataProducts/inc/Doublet.hh"

//ROOT
#include "TStopwatch.h"

namespace fhicl {
  class ParameterSet;
}

// class mu2e::Tracker;
// class mu2e::TrackerCalibrations;

namespace mu2e {
  class KalFitHack : public KalContext {
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
    int                         _debug;
    unsigned                    _minnstraws;
    double                      _maxmatfltdiff; // maximum difference in track flightlength to separate to intersections of the same material
    vector<bool>                _weedhits;
    double                      _maxhitchi;
    unsigned                    _maxweed;
    std::vector<double>         _hiterr;
    double                      _maxdriftpull;
    fhicl::ParameterSet*        _darPset;         // parameter set for doublet ambig resolver
    std::vector<AmbigResolver*> _ambigresolver;
    bool                        _initt0;
    bool                        _updateT0;
    int                         _updateT0Mode;    // 0: use cluster T0 1:update T0 assuming no cluster time
    double                      _minHitDrift;
    double                      fRdriftMinusDocaTol;
    int                         fSign[4][2];
    int                         _daveMode;
    std::vector<double>         _t0tol;
    TStopwatch*                 fStopwatch;       // = new TStopwatch();
    double                      _t0errfac;       // fudge factor for the calculated t0 error
    double                      _mint0doca;      // minimum (?) doca for t0 hits
    double                      _t0nsig;	        // # of sigma to include when selecting hits for t0
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

    const mu2e::PtrStepPointMCVectorCollection*  fListOfMCStrawHits;

    const mu2e::Tracker*             _tracker;     // straw tracker geometry
    const mu2e::TrackerCalibrations* _tcal;
//-----------------------------------------------------------------------------
// constructors and destructor, parameter set should be passed in on construction
//-----------------------------------------------------------------------------
  public:
    explicit KalFitHack(fhicl::ParameterSet const&);
    virtual ~KalFitHack();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
    int  nIter       () { return _nIter        ; } 
    int  maxIteration() { return _hiterr.size()-1; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
    void setNIter        (int N   ) { _nIter        = N  ; }

    void setTracker      (const Tracker*             Tracker) { _tracker = Tracker; }
    void setTrackerCalib (const TrackerCalibrations* TCal   ) { _tcal    = TCal;    }

    void setStepPointMCVectorCollection(const mu2e::PtrStepPointMCVectorCollection* List) {
      fListOfMCStrawHits = List;
    }
//-----------------------------------------------------------------------------
// add a set of hits to an existing fit
//-----------------------------------------------------------------------------
    virtual void addHits(KalFitResult&             kres   , 
			 const StrawHitCollection* straws , 
			 std::vector<hitIndex>     indices, 
			 double                    maxchi ,
			 CalTimePeak*              TPeak=NULL );

    void findBoundingHits(std::vector<TrkStrawHit*>&                   hits, 
			  double                                       flt0,
			  std::vector<TrkStrawHit*>::reverse_iterator& ilow ,
			  std::vector<TrkStrawHit*>::iterator&         ihigh);

    bool fitable     (TrkDefHack const& tdef);
    void fitIteration(KalFitResult& kres , int Iteration, CalTimePeak* TPeak);
    void fitTrack    (KalFitResult& kres , CalTimePeak* TPeak=NULL);
    void initCaloT0  (CalTimePeak*  TPeak, TrkDefHack const& tdef, TrkT0& t0);
    void initT0      (TrkDefHack const& tdef , TrkT0& t0);
//-----------------------------------------------------------------------------
// main function: given a track definition, create a fit object from it
//-----------------------------------------------------------------------------
    virtual void makeTrack      (KalFitResult& kRes, CalTimePeak* TPeak=NULL);
//---------------------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following function for printing the hits included in the track
//----------------------------------------------------------------------------------------------
    void printHits       (KalFitResult& kres, const char* Caller);

    bool unweedHits      (KalFitResult& kres, double maxchi);
    void updateCalT0     (KalFitResult& kres, CalTimePeak* TPeak);
    void updateHitTimes  (KalFitResult& kres);
    bool updateT0        (KalFitResult& kres);
    bool weedHits        (KalFitResult& kres, int Iteration);

// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
//-----------------------------------------------------------------------------
// overloaded functions of KalContext
//-----------------------------------------------------------------------------
    virtual void makeHits     (KalFitResult& kres, TrkT0 const& t0);
    virtual void makeMaterials(KalFitResult& kres);
//-----------------------------------------------------------------------------
// latest code by Dave
//-----------------------------------------------------------------------------
    double   zFlight    (KalRep* KRep, double pz); // *** FIXME category
    unsigned addMaterial(KalRep* KRep);

  };
}
#endif
