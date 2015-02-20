//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitHack.hh,v 1.3 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
//
#ifndef KalFitHack_HH
#define KalFitHack_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
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

namespace mu2e {
  class KalFitHack : public KalContext {
  public:
// define different ambiguity resolution strategies
    enum ambigStrategy {fixedambig=0,pocaambig=1,hitambig=2,panelambig=3,hitambig_hack=4};
// parameter set should be passed in on construction
    explicit KalFitHack(fhicl::ParameterSet const&);
    virtual ~KalFitHack();
//-----------------------------------------------------------------------------
// main function: given a track definition, create a fit object from it
//-----------------------------------------------------------------------------
    virtual void makeTrack(KalFitResult& kres, CalTimePeak* TPeak=NULL, int markDoubs=0);

//---------------------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following function for printing the hits included in the track
//----------------------------------------------------------------------------------------------
    void printHits(KalFitResult& kres);

//----------------------------------------------------------------------    
// 2015 - 02 - 17 G. Pezzullo added a method for searching the doublets
    void findDoublets      (KalRep* krep,   std::vector<TrkStrawHit*> *hits,  DoubletCollection *dcol);//search doublets in a giventimepeak
    void markDoublets      (DoubletCollection *dcol);
    void findAndMarkMultiplets(KalRep* krep,   std::vector<TrkStrawHit*> *hits);
    
    void setMultipletsAmbig(std::vector<TrkStrawHit*> *hits);

//-----------------------------------------------------------------------------------------
//2015 - 02 -20 G. Pezzu addedd the function for calculataing the slope of the lines
// tangent to two given circles
//-----------------------------------------------------------------------------------------
    void findLines(double xa, double ya, double ra,
		   double xb, double yb, double rb,
		   double *slopes);
    

// add a set of hits to an existing fit
    virtual void addHits(KalFitResult&             kres   , 
			 const StrawHitCollection* straws , 
			 std::vector<hitIndex>     indices, 
			 double                    maxchi );

// Try to put back inactive hits
    bool unweedHits(KalFitResult& kres, double maxchi);
    int  NIter                () { return fNIter; }
    void SetNIter             (int N) { fNIter = N; }
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
  protected:
    // configuration parameters
    int                         _debug;
    bool                        _weedhits;
    double                      _maxhitchi;
    unsigned                    _maxweed;
    std::vector<double>         _herr;
    double                      _maxdriftpull;
    std::vector<AmbigResolver*> _ambigresolver;
    bool                        _initt0;
    bool                        _updatet0;
    double                      fMinHitDrift;
    double                      fRdriftMinusDocaTol;
    int                         _daveMode;
    std::vector<double>         _t0tol;

    bool fitable              (TrkDef const& tdef);
    bool weedHits             (KalFitResult& kres);
    void initCaloT0           (CalTimePeak* TPeak, TrkDef const& tdef, TrkT0& t0);
    void initT0               (TrkDef const& tdef, TrkT0& t0);
    void updateCalT0          (KalFitResult& kres, CalTimePeak* TPeak);
    bool updateT0             (KalFitResult& kres);
    void fitTrack             (KalFitResult& kres, CalTimePeak* TPeak=NULL);
//-----------------------------------------------------------------------------
// overloaded functions of KalContext
//-----------------------------------------------------------------------------
    virtual void makeHits     (KalFitResult& kres, TrkT0 const& t0);
    virtual void makeMaterials(KalFitResult& kres);

  private:
    TStopwatch*        fStopwatch;       // = new TStopwatch();
    double             _t0errfac;       // fudge factor for the calculated t0 error
    double             _mint0doca;      // minimum (?) doca for t0 hits
    double             _t0nsig;	        // # of sigma to include when selecting hits for t0
    double             _dtoffset;       // track - luster time offset, ns
    double             fScaleErrDoublet;
    int                fMarkDoublets;
    int                fILoopMarkDoublets;
    bool               _removefailed;
    unsigned           _minnstraws;
    TrkParticle        _tpart;
    TrkFitDirection    _fdir;
    std::vector<int>   _ambigstrategy;
    mutable BField*    _bfield;
    int                 fNIter;
    const CalTimePeak*  fTimePeak;
    int                 fAmbigVec[40000];
    int                 fAmbigVecSlope[40000];
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    void fitIteration    (KalFitResult& kres, size_t iiter, CalTimePeak* TPeak=NULL);

    void updateHitTimes  (KalFitResult& kres);

    void findBoundingHits(std::vector<TrkStrawHit*>&                   hits, 
			  double                                       flt0,
			  std::vector<TrkStrawHit*>::reverse_iterator& ilow ,
			  std::vector<TrkStrawHit*>::iterator&         ihigh);
  };
}
#endif
