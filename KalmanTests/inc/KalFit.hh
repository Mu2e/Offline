//
// Object to perform BaBar Kalman fit
//
// $Id: KalFit.hh,v 1.21 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/23 17:52:27 $
//
#ifndef KalFit_HH
#define KalFit_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
// KalFit objects
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/AmbigResolver.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "BField/BField.hh"
#include "TrkBase/TrkParticle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++

namespace mu2e 
{
// struct defining the output of the fit
  struct TrkKalFit {
    TrkT0 _t00; // t0 value at start of fit
    TrkT0 _t0; // t0 value at end of fit
    KalRep* _krep; // Kalman rep, owned by the trk
    std::vector<TrkStrawHit*> _hits; // straw hits, owned by the KalRep
    TrkErrCode _fit; // error code from last fit
    unsigned _nt0iter; // number of times t0 was iterated
    unsigned _nweediter; // number of iterations on hit weeding
    TrkKalFit() : _krep(0), _fit(TrkErrCode::fail) {}
    ~TrkKalFit() { delete _krep;}
    void setT0(const TrkT0& t0) { _t0 = t0; }
    void setT00(const TrkT0& t0) { _t00 = t0; }
    void removeFailed() { if(_fit.failure())deleteTrack(); }
    void fit() { if(_fit.success()) _fit = _krep->fit(); }
    void deleteTrack();
    KalRep* stealTrack() { KalRep* retval = _krep; _krep=0; _hits.clear(); return retval; }
  };
  
  class KalFit
  {
  public:
// define different t0 strategies.  Eventually t0 finding should be its own class
    enum t0Strategy {external=0,median,histogram};
// define different ambiguity resolution strategies.  These will eventually be their own classes
    enum ambigStrategy {fixedambig=0,pocaambig,hitambig,panelambig};
// parameter set should be passed in on construction
    explicit KalFit(fhicl::ParameterSet const&);
    virtual ~KalFit();
// main function: given a track definition, create a fit object from it
    void makeTrack(TrkDef const& mytrk,TrkKalFit& myfit);
// add a set of hits to an existing fit
    void addHits(TrkKalFit& myfit,const StrawHitCollection* straws, std::vector<hitIndex> indices);
  private:
// Fetch the BField.  this function fetches the field if it's not yet initialized
    const BField* bField();
    KalContext* _kalcon; // BaBar configuration object
    BField* _bfield; // magentic field description, BaBar wrapper, OWNED BY THIS CLASS
    // configuration parameters
    int _debug;
    bool _fieldcorr;
    bool _material;
//    bool _ambigflip;
    bool _weedhits;
    bool _updatet0;
    bool _removefailed;
    double _t0tol;
    double _maxhitchi;
    unsigned _maxiter;
    double _mingap;
    unsigned _minnstraws;
    unsigned _minndof;
    unsigned _maxweed;
    std::vector<double> _herr;
    double _ssmear;
    double _t0errfac;
    double _mint0doca;
    double _maxdriftpull;
    double _t0nsig;
    TrkParticle _tpart;
    t0Strategy _t0strategy;
    std::vector<int> _ambigstrategy;
    std::vector<AmbigResolver*> _ambigresolver;
    // helper functions
    bool fitable(TrkDef const& mytrk);
    bool updateT0(TrkKalFit& myfit);
    bool weedHits(TrkKalFit& myfit);
    void fitTrack(TrkKalFit& myfit);
    void makeHits(TrkDef const& mytrk,TrkKalFit& myfit);
    void initT0(TrkDef const& mytrk,TrkT0& t0);
    double findZFltlen(const TrkKalFit& myfit,double zval);
// general
    static const double _vlight;
// 
  };
}
#endif
