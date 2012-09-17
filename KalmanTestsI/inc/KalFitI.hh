//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitI.hh,v 1.2 2012/09/17 14:44:30 ignatov Exp $
// $Author: ignatov $ 
// $Date: 2012/09/17 14:44:30 $
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
#include "KalmanTests/inc/KalFitResult.hh"
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
    void makeTrack(KalFitResult& kdef);
// add a set of hits to an existing fit
    void addHits(KalFitResult& kdef,std::vector<hitIndex> indices,bool active=true);
    void addHitsUnique(KalFitResult& kdef,std::vector<hitIndex> indices,bool active=true);

    void reActivateHitsbyTurn(KalFitResult& kdef);
    void reActivateHitsbyChi2(KalFitResult& kdef);

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
    bool updateT0(KalFitResult& kdef);
    bool weedHits(KalFitResult& kdef);
    bool unweedHits(KalFitResult& kdef);
    void fitTrack(KalFitResult& kdef);
    void makeHits(KalFitResult& kdef,TrkT0 const& t0);
    bool fixHitTurn(KalFitResult& kres,TrkStrawHit* trkhit);
    void initT0(TrkDef const& mytrk,TrkT0& t0);
    double findZFltlen(KalFitResult& kdef,double zval);
// general
    static const double _vlight;
    static const double _vdrift;
// 
  };
}
#endif
