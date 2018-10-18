//
// Track definition object
//
// $Id: TrkDefHack.cc,v 1.9 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
// Original author David Brown, LBNL
//
#include "BTrk/BaBar/BaBar.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/TrkDefHack.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;

namespace mu2e 
{
    // static variables
  HepVector       TrkDefHack::_dpar(5,0);
  HepSymMatrix    TrkDefHack::_dcov(5,0);
  TrkParticle     TrkDefHack::_eminus(TrkParticle::e_minus);
  TrkFitDirection TrkDefHack::_downstream(TrkFitDirection::downstream);


//-----------------------------------------------------------------------------
// ignore helix part, initialize only hit lists
//-----------------------------------------------------------------------------
  TrkDefHack::TrkDefHack(const KalSeed*                     Seed,
			 const StrawHitCollection*          StrawCollection ,
			 const StrawHitPositionCollection*  ShPosCollection , 
			 const StrawHitFlagCollection*      ShFlagCollection,
			 TrkParticle     const&             tpart,
			 TrkFitDirection const&             fdir ):
    _h0(_dpar,_dcov),
    _tpart(tpart),
    _fdir(fdir)
  {
    int nhits = Seed->hits().size();
    for (int i=0; i<nhits; ++i) {
      const mu2e::TrkStrawHitSeed*  hit = &Seed->hits().at(i);
      if (!hit->_flag.hasAnyProperty(StrawHitFlag::outlier)){
	 _timeCluster._strawHitIdxs.push_back(hit->index());
      }
    }

    _shcol  = StrawCollection;
    _shpos  = ShPosCollection;
    _shfcol = ShFlagCollection;
  }

//-----------------------------------------------------------------------------
// ignore helix part, initialize only hit lists
//-----------------------------------------------------------------------------
  TrkDefHack::TrkDefHack(const HelixSeed*                   Seed,
			 const StrawHitCollection*          StrawCollection ,
			 const StrawHitPositionCollection*  ShPosCollection , 
			 const StrawHitFlagCollection*      ShFlagCollection,
			 TrkParticle     const&             tpart,
			 TrkFitDirection const&             fdir ):
    _h0(_dpar,_dcov),
    _tpart(tpart),
    _fdir(fdir)
  {

    // get the collection of the StrawHitIndices

    int nHits = Seed->hits().size();
    for (int i=0; i<nHits; ++i){
      const HelixHit* hit = &Seed->hits().at(i);
      if (!hit->_flag.hasAnyProperty(StrawHitFlag::outlier)){
	 _timeCluster._strawHitIdxs.push_back(hit->index());
      }
    }

    _shcol  = StrawCollection;
    _shpos  = ShPosCollection;
    _shfcol = ShFlagCollection;
  }

//-----------------------------------------------------------------------------
  TrkDefHack::TrkDefHack(const CalTimePeak*                 TimePeak,
			 const StrawHitCollection*          StrawCollection ,
			 const StrawHitPositionCollection*  ShPosCollection , 
			 const StrawHitFlagCollection*      ShFlagCollection,
			 TrkParticle     const&             tpart,
			 TrkFitDirection const&             fdir ):
    _h0(_dpar,_dcov),
    _tpart(tpart),
    _fdir(fdir)
  {
    _timeCluster._strawHitIdxs = TimePeak->_index;

    _shcol  = StrawCollection;
    _shpos  = ShPosCollection;
    _shfcol = ShFlagCollection;
  }

//-----------------------------------------------------------------------------
  TrkDefHack::TrkDefHack(const TimeCluster*                 TimeCluster,
			 const StrawHitCollection*          StrawCollection ,
			 const StrawHitPositionCollection*  ShPosCollection , 
			 const StrawHitFlagCollection*      ShFlagCollection,
			 TrkParticle     const&             tpart,
			 TrkFitDirection const&             fdir ):
    _h0(_dpar,_dcov),
    _tpart(tpart),
    _fdir(fdir)
  {
    if (TimeCluster) _timeCluster = *TimeCluster;

    _shcol  = StrawCollection;
    _shpos  = ShPosCollection;
    _shfcol = ShFlagCollection;
  }

//-----------------------------------------------------------------------------
// shouldn't 'Hits' be coming as a part of 'TimeCluster' ?
//-----------------------------------------------------------------------------
  TrkDefHack::TrkDefHack(const HelixTraj&                   Helix,
			 const TimeCluster*                 TimeCluster,
			 const StrawHitCollection*          StrawCollection ,
			 const StrawHitPositionCollection*  ShPosCollection , 
			 const StrawHitFlagCollection*      ShFlagCollection,
			 TrkParticle     const&             tpart,
			 TrkFitDirection const&             fdir ):
    _h0(Helix),
    _tpart(tpart),
    _fdir(fdir)
  {
    if (TimeCluster) _timeCluster = *TimeCluster;
    
    _shcol  = StrawCollection;
    _shpos  = ShPosCollection;
    _shfcol = ShFlagCollection;
  }

//-----------------------------------------------------------------------------
  TrkDefHack::~TrkDefHack(){
  }
//-----------------------------------------------------------------------------
// invalidate pointers to hits and flags
//-----------------------------------------------------------------------------
  void TrkDefHack::init() {
    _shcol  = NULL;
    _shpos  = NULL;
    _shfcol = NULL;
    _timeCluster._strawHitIdxs.clear();
  }
}
