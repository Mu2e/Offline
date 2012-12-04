//
// Object to perform BaBar Kalman fit
//
// $Id: KalFit.hh,v 1.27 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/12/04 00:51:26 $
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
// C++

namespace mu2e 
{
  class KalFit : public KalContext
  {
  public:
// define different ambiguity resolution strategies
    enum ambigStrategy {fixedambig=0,pocaambig,hitambig,panelambig};
// parameter set should be passed in on construction
    explicit KalFit(fhicl::ParameterSet const&);
    virtual ~KalFit();
// main function: given a track definition, create a fit object from it
    virtual void makeTrack(KalFitResult& kres);
// add a set of hits to an existing fit
    virtual void addHits(KalFitResult& kres,const StrawHitCollection* straws, std::vector<hitIndex> indices, double maxchi);
// Try to put back inactive hits
    bool unweedHits(KalFitResult& kres, double maxchi);
// KalContext interface
    virtual const TrkVolume* trkVolume(trkDirection trkdir) const ;
    BField const& bField() const;
  protected:
    // configuration parameters
    int _debug;
    bool _weedhits;
    double _maxhitchi;
    unsigned _maxweed;
    std::vector<double> _herr;
    double _maxdriftpull;
    std::vector<AmbigResolver*> _ambigresolver;
    bool _initt0;
    bool _updatet0;
    std::vector<double> _t0tol;
    bool fitable(TrkDef const& tdef);
    bool weedHits(KalFitResult& kres);
    void initT0(TrkDef const& tdef, TrkT0& t0);
    bool updateT0(KalFitResult& kres);

  private:
    double _t0errfac; // fudge factor for the calculated t0 error
    double _mint0doca; // minimum (?) doca for t0 hits
    double _t0nsig; // # of sigma to include when selecting hits for t0
    bool _removefailed;
    unsigned _minnstraws;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    std::vector<int> _ambigstrategy;
    mutable BField* _bfield;
    // helper functions

    virtual void fitTrack(KalFitResult& kres);
    void fitIteration(KalFitResult& kres,size_t iiter);
    virtual void makeHits(KalFitResult& kres, TrkT0 const& t0);
    virtual void makeMaterials(KalFitResult& kres);
    void updateHitTimes(KalFitResult& kres);
    void findBoundingHits(std::vector<TrkStrawHit*>& hits, double flt0,
	std::vector<TrkStrawHit*>::reverse_iterator& ilow,
	std::vector<TrkStrawHit*>::iterator& ihigh);
  };
}
#endif
