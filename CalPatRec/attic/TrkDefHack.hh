//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDefHack.hh,v 1.18 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/30 12:19:38 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDefHack_HH
#define TrkDefHack_HH
// Mu2e
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
// BTrk includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkT0.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

namespace mu2e  {
  
  class CalTimePeak;
  class HelixSeed;
  class KalSeed;
  
  class TrkDefHack {
  public:
    HelixTraj                  _h0;          // helix estimate, valid in the region around z=0
    TimeCluster                _timeCluster; // t0 and hit indices
    TrkParticle                _tpart;       // particle type. Defines both the charge and the mass
    TrkFitDirection            _fdir;        // fit direction

    const StrawHitCollection*         _shcol;
    const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _shfcol;

    // dummy variables

    static CLHEP::HepVector    _dpar;
    static CLHEP::HepSymMatrix _dcov;
    static TrkParticle         _eminus;
    static TrkFitDirection     _downstream;

  public:

    TrkDefHack(const KalSeed*                     Seed                   ,
	       const StrawHitCollection*          StrawCollection  = NULL,
	       const StrawHitPositionCollection*  ShPosCollection  = NULL, 
	       const StrawHitFlagCollection*      ShFlagCollection = NULL, 
	       TrkParticle     const&             tpart = _eminus    ,
	       TrkFitDirection const&             fdir  = _downstream);

    TrkDefHack(const HelixSeed*                   Seed                   ,
	       const StrawHitCollection*          StrawCollection  = NULL,
	       const StrawHitPositionCollection*  ShPosCollection  = NULL, 
	       const StrawHitFlagCollection*      ShFlagCollection = NULL, 
	       TrkParticle     const&             tpart = _eminus    ,
	       TrkFitDirection const&             fdir  = _downstream);

    TrkDefHack(const CalTimePeak*                 TimePeak               ,
	       const StrawHitCollection*          StrawCollection  = NULL,
	       const StrawHitPositionCollection*  ShPosCollection  = NULL, 
	       const StrawHitFlagCollection*      ShFlagCollection = NULL, 
	       TrkParticle     const&             tpart = _eminus    ,
	       TrkFitDirection const&             fdir  = _downstream);

    TrkDefHack(const TimeCluster*                 TCluster         = NULL,
	       const StrawHitCollection*          StrawCollection  = NULL,
	       const StrawHitPositionCollection*  ShPosCollection  = NULL, 
	       const StrawHitFlagCollection*      ShFlagCollection = NULL, 
	       TrkParticle     const&             tpart = _eminus    ,
	       TrkFitDirection const&             fdir  = _downstream);

    TrkDefHack(HelixTraj       const&             Helix	      ,
	       const TimeCluster*                 TCluster         = NULL,
	       const StrawHitCollection*          StrawCollection  = NULL,
	       const StrawHitPositionCollection*  ShPosCollection  = NULL, 
	       const StrawHitFlagCollection*      ShFlagCollection = NULL, 
	       TrkParticle     const&             tpart = _eminus    ,
	       TrkFitDirection const&             fdir  = _downstream);

    ~TrkDefHack();

    // accessors

    const std::vector<StrawHitIndex>& strawHitIndices() const { return _timeCluster._strawHitIdxs;}
    const HelixTraj&            helix       () const { return _h0; }
    const CLHEP::HepSymMatrix&  helixCovMatr() const { return _h0.parameters()->covariance(); }
    const CLHEP::HepSymMatrix&  covMatr     () const { return _dcov; }
    TrkParticle const&          particle    () const { return _tpart; }
    TrkFitDirection const&      fitdir      () const { return _fdir; }
    TrkT0 const&                t0          () const { return _timeCluster._t0; }
    
    void setT0     (TrkT0     const& t0   ) { _timeCluster._t0 = t0; }
    void setHelix  (HelixTraj const& helix) { _h0 = helix; }

    // this is a vector copy by value - why ?  
    void setIndices(std::vector<StrawHitIndex> const& indices ) { _timeCluster._strawHitIdxs = indices; }
    int  hitIndex(int I) { return _timeCluster._strawHitIdxs[I]; }

    // append a straw hit to this track definition
    void appendHit(size_t index) { _timeCluster._strawHitIdxs.push_back(index); }
    
    const StrawHitPositionCollection* shpos () const { return _shpos ; }
    const StrawHitFlagCollection*     shfcol() const { return _shfcol; }
    const StrawHitCollection*         shcol () const { return _shcol;  }

    void init();
  };
}

#endif
