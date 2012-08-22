//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitI.hh,v 1.1 2012/08/22 17:30:37 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/08/22 17:30:37 $
//
#ifndef KalFitI_HH
#define KalFitI_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
//#include "BaBar/PdtPid.hh"
//#include "TrkBase/TrkRecoTrk.hh"
// KalFitI objects
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/AmbigResolver.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "BField/BField.hh"
#include "KalmanTests/inc/KalFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++
#include <iostream>
#include <vector>

#include "TrkBase/TrkParticle.hh"

namespace mu2e 
{
  class KalFitI
  {
  public:
// define the fit direction as downstream (towards positive Z) or upstream (towards negative Z).
    enum fitDirection {downstream=0,upstream};
// define different t0 strategies.  Eventually t0 finding should be its own class
    enum t0Strategy {external=0,median,histogram};
// parameter set should be passed in on construction
    explicit KalFitI(fhicl::ParameterSet const&);
    virtual ~KalFitI();
// main function: given a track definition, create a fit object from it
    void makeTrack(TrkDef& mytrk,TrkKalFit& myfit);
// add a set of hits to an existing fit
    void addHits(TrkDef& mytrk,TrkKalFit& myfit,std::vector<hitIndex> indices,bool active=true);
    void addHitsUnique(TrkDef& mytrk,TrkKalFit& myfit,std::vector<hitIndex> indices,bool active=true);

    void reActivateHitsbyTurn(TrkDef& mytrk,TrkKalFit& myfit);
    void reActivateHitsbyChi2(TrkDef& mytrk,TrkKalFit& myfit);

// sim helper
    double _flt0;
    std::vector<double> _hitflt;

    const BField* bField();
  private:
// Fetch the BField.  this function fetches the field if it's not yet initialized
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
//    double _herr;
    std::vector<double> _herr;
    double _ssmear;
    double _t0errfac;
    double _mint0doca;
    double _maxdriftpull;
    double _t0nsig;
    TrkParticle _tpart;
    t0Strategy _t0strategy;
    bool _exacthitturn;
    std::vector<int> _ambigstrategy;
    std::vector<AmbigResolver*> _ambigresolver;
    // helper functions
    bool fitable(TrkDef const& mytrk);
    bool updateT0(TrkKalFit& myfit);
    bool weedHits(TrkKalFit& myfit);
    bool unweedHits(TrkKalFit& myfit);
    void fitTrack(TrkKalFit& myfit);
    void makeHits(TrkDef const& mytrk,TrkKalFit& myfit);
    bool fixHitTurn(const TrkDef& mytrk,TrkStrawHit* trkhit);
    void initT0(TrkDef const& mytrk,TrkT0& t0);
    double findZFltlen(const TrkKalFit& myfit,double zval);
// general
    static const double _vlight;
    static const double _vdrift;
// 
  };
}
#endif
