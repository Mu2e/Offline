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
#include "BaBar/BaBar.hh"
// KalFitHack objects
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/AmbigResolver.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "BField/BField.hh"
#include "TrkBase/TrkParticle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/Doublet.hh"

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
    std::vector<AmbigResolver*> _ambigresolver;
    bool                        _initt0;
    bool                        _updatet0;
    double                      fMinHitDrift;
    double                      fRdriftMinusDocaTol;
    int                         fSign[4][2];
    int                         _daveMode;
    std::vector<double>         _t0tol;
    std::vector<Doublet>        fListOfDoublets;
    TStopwatch*                 fStopwatch;       // = new TStopwatch();
    double                      _t0errfac;       // fudge factor for the calculated t0 error
    double                      _mint0doca;      // minimum (?) doca for t0 hits
    double                      _t0nsig;	        // # of sigma to include when selecting hits for t0
    double                      _dtoffset;       // track - luster time offset, ns
    double                      fScaleErrDoublet;
    int                         fUseDoublets;
    int                         fILoopUseDoublets;
    double                      fMinDriftDoublet;
    double                      fDeltaDriftDoublet;
    double                      fSigmaSlope;
    std::string                 fMakeStrawHitModuleLabel;
		                
    bool                        _removefailed;
    unsigned                    _minnstraws;
    TrkParticle                 _tpart;
    TrkFitDirection             _fdir;
    std::vector<int>            _ambigstrategy;
    mutable BField*             _bfield;
    int                         fNIter;
    const CalTimePeak*          fTimePeak;
    int                         fAmbigVec     [40000];
    int                         fAmbigVecSlope[40000];
    int                         fAnnealingStep;
    int                         fDecisionMode; // 0:decision is not forced; 1:decision has to be made

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
    std::vector<Doublet>*  listOfDoublets() { return &fListOfDoublets; }
    int                    NIter         () { return fNIter;           }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
    void setDecisionMode (int Mode) { fDecisionMode = Mode; }
    void SetNIter        (int N   ) { fNIter        = N  ; }

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
//----------------------------------------------------------------------    
// 2015-02-17 G.Pezzullo: search doublets in a given timepeak
//-----------------------------------------------------------------------------
    void          findDoublets(KalFitResult& KRes, DoubletCollection *DCol);
//-----------------------------------------------------------------------------------------
// 2015-02-20: G.Pezzu added function for calculataing the slope of the lines
// tangent to two given circles
//-----------------------------------------------------------------------------------------
    void findLines(Hep3Vector A[2], double rb[2], double *Slopes);

    const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& Hits, 
				const KalRep*                    Krep) const  ;

    bool fitable     (TrkDef const& tdef);
    void fitIteration(KalFitResult& kres , size_t iiter, CalTimePeak* TPeak=NULL);
    void fitTrack    (KalFitResult& kres , CalTimePeak* TPeak=NULL);
    void initCaloT0  (CalTimePeak*  TPeak, TrkDef const& tdef, TrkT0& t0);
    void initT0      (TrkDef const& tdef , TrkT0& t0);
//-----------------------------------------------------------------------------
// main function: given a track definition, create a fit object from it
//-----------------------------------------------------------------------------
    virtual void makeTrack      (KalFitResult& kres, CalTimePeak* TPeak=NULL, int markDoubs=0);
    void         markMultiplets (KalFitResult& Kres, DoubletCollection *dcol);
    void         markDoublet    (Doublet *doub, int index0, int index1);
//---------------------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following function for printing the hits included in the track
//----------------------------------------------------------------------------------------------
    void printHits(KalFitResult& kres, const char* Caller);

    bool unweedHits    (KalFitResult& kres, double maxchi);
    void updateCalT0   (KalFitResult& kres, CalTimePeak* TPeak);
    void updateHitTimes(KalFitResult& kres);
    bool updateT0      (KalFitResult& kres);
    bool weedHits      (KalFitResult& kres);

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
