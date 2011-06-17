//
// Object to perform BaBar Kalman fit
//
// $Id: KalFit.hh,v 1.4 2011/06/17 21:56:57 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/17 21:56:57 $
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
#include "KalmanTests/inc/DetStrawHitElem.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "BField/BField.hh"
#include "MatEnv/MatDBInfo.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++

using namespace std; 

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
    TrkKalFit() : _trk(0),_krep(0) {}
    ~TrkKalFit() { delete _trk;}
    void setT0(const TrkT0& t0) { _t0 = t0; }
    void setT00(const TrkT0& t0) { _t00 = t0; }
    void removeFailed() { if(_fit.failure())deleteTrack(); }
    void fit() { _fit = _krep->fit(); }
    void deleteTrack() { delete _trk; _trk=0; _krep = 0; }
  };
  
  class KalFit
  {
  public:
// parameter set should be passed in on construction
    explicit KalFit(fhicl::ParameterSet const&);
    virtual ~KalFit();
// main function: given a track definition, create a fit object from it
    void makeTrack(TrkDef const& mytrk,TrkKalFit& myfit);
  private:
// Fetch the BField.  this function fetches the field if it's not yet initialized
    const BField* bField();
    KalContext* _kalcon; // BaBar configuration object
    BField* _bfield; // magentic field description, BaBar wrapper, OWNED BY THIS CLASS
    DetStrawHitElem _wallelem; // fake element to represent straw hit material
    DetStrawHitElem _gaselem; // fake element to represent straw hit material
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
    unsigned _minnstraws;
    unsigned _maxweed;
    // helper functions
    bool fitable(TrkDef const& mytrk);
    void makeHits(TrkDef const& mytrk, TrkKalFit& myfit);
    bool updateT0(TrkKalFit& myfit);
    bool weedHits(TrkKalFit& myfit);
// general
    static const double _vlight;
    static const double _vdrift;
    static MatDBInfo* _matdbinfo;
// 
  };
}
#endif
