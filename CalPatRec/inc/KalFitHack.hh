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
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
// KalFitHack objects
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/AmbigResolver.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "KalmanTests/inc/Doublet.hh"

//ROOT
#include "TStopwatch.h"

namespace fhicl {
  class ParameterSet;
}

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
    bool                        _weedhits;
    double                      _maxhitchi;
    unsigned                    _maxweed;
    std::vector<double>         _hiterr;
    double                      _maxdriftpull;
    fhicl::ParameterSet*        _darPset;         // parameter set for doublet ambig resolver
    std::vector<AmbigResolver*> _ambigresolver;
    bool                        _initt0;
    bool                        _updateT0;
    int                         _updateT0Mode;    // 0: use cluster T0 1:update T0 assuming no cluster time
    double                      fMinHitDrift;
    double                      fRdriftMinusDocaTol;
    int                         fSign[4][2];
    int                         _daveMode;
    std::vector<double>         _t0tol;
    //    std::vector<Doublet>        fListOfDoublets;
    TStopwatch*                 fStopwatch;       // = new TStopwatch();
    double                      _t0errfac;       // fudge factor for the calculated t0 error
    double                      _mint0doca;      // minimum (?) doca for t0 hits
    double                      _t0nsig;	        // # of sigma to include when selecting hits for t0
    double                      _dtoffset;       // track - luster time offset, ns
    double                      fScaleErrDoublet;
    //    int                         fUseDoublets;
    double                      fMinDriftDoublet;
    double                      fDeltaDriftDoublet;
    double                      _maxDoubletChi2;
    double                      fSigmaSlope;
    std::string                 fMakeStrawHitModuleLabel;
		                
    bool                        _removefailed;
    unsigned                    _minnstraws;
    TrkParticle                 _tpart;
    TrkFitDirection             _fdir;
    std::vector<int>            _ambigstrategy;
    mutable BField*             _bfield;
    int                         _nIter;
    const CalTimePeak*          fTimePeak;
    int                         _annealingStep;

    const mu2e::PtrStepPointMCVectorCollection*  fListOfMCStrawHits;
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

//     const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& Hits, 
// 				const KalRep*                    Krep) const  ;

    bool fitable     (TrkDef const& tdef);
    void fitIteration(KalFitResult& kres , int Iteration, CalTimePeak* TPeak, int Final);
    void fitTrack    (KalFitResult& kres , CalTimePeak* TPeak=NULL);
    void initCaloT0  (CalTimePeak*  TPeak, TrkDef const& tdef, TrkT0& t0);
    void initT0      (TrkDef const& tdef , TrkT0& t0);
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
    bool weedHits        (KalFitResult& kres, int Iteration, int Final);

// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
//-----------------------------------------------------------------------------
// overloaded functions of KalContext
//-----------------------------------------------------------------------------
    virtual void makeHits     (KalFitResult& kres, TrkT0 const& t0);
    virtual void makeMaterials(KalFitResult& kres);
  };
}
#endif
