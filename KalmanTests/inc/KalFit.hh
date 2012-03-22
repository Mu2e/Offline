//
// Object to perform BaBar Kalman fit
//
// $Id: KalFit.hh,v 1.16 2012/03/22 22:32:23 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/03/22 22:32:23 $
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
#include "BaBar/PdtPid.hh"
#include "TrkBase/TrkRecoTrk.hh"
// KalFit objects
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "BField/BField.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++

namespace mu2e 
{
// struct defining the output of the fit
  struct TrkKalFit {
    TrkT0 _t00; // t0 value at start of fit
    TrkT0 _t0; // t0 value at end of fit
    TrkRecoTrk* _trk; // the track
    KalRep* _krep; // Kalman rep, owned by the trk
    std::vector<TrkStrawHit*> _hits; // straw hits, owned by the KalRep
    TrkErrCode _fit; // error code from last fit
    unsigned _nt0iter; // number of times t0 was iterated
    unsigned _nweediter; // number of iterations on hit weeding
    TrkKalFit() : _trk(0),_krep(0), _fit(TrkErrCode::fail) {}
    ~TrkKalFit() { delete _trk;}
    void setT0(const TrkT0& t0) { _t0 = t0; }
    void setT00(const TrkT0& t0) { _t00 = t0; }
    void removeFailed() { if(_fit.failure())deleteTrack(); }
    void fit() { _fit = _krep->fit(); }
    void deleteTrack();
    TrkRecoTrk* stealTrack() { TrkRecoTrk* retval = _trk; _trk=0; _krep=0; _hits.clear(); return retval; }
  };
  
  class KalFit
  {
  public:
// define different t0 strategies.  Eventually t0 finding should be its own class
    enum t0Strategy {external=0,median,histogram};
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
    bool _ambigflip;
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
    double _herr;
    double _ssmear;
    double _t0errfac;
    double _mint0doca;
    double _maxdriftpull;
    double _t0nsig;
    int _fitpart;
    t0Strategy _t0strategy;
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
    static const double _vdrift;
// 
  };
}
#endif
