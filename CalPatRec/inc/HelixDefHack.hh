//
// Object to perform helix fit to straw hits
//
// $Id: HelixDefHack.hh,v 1.2 2014/05/18 13:56:50 murat Exp $
// $Author: murat $ 
// $Date: 2014/05/18 13:56:50 $
//
#ifndef HelixDefHack_HH
#define HelixDefHack_HH

// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

// #include "RecoDataProducts/inc/StrawHitFlag.hh"

#include "KalmanTests/inc/TrkDef.hh"

namespace mu2e {
					// add the StrawHitPosition collection to TrkDef
  class HelixDefHack : public TrkDef {
  private:
    const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _shfcol;
  public:
    HelixDefHack(TrkDef const& tdef) : TrkDef(tdef), _shpos(0) {}

    HelixDefHack(const StrawHitCollection*         strawcollection,
		 const StrawHitPositionCollection* shposcollection, 
		 const StrawHitFlagCollection*     ShFlagCollection, 
		 const std::vector<hitIndex>&      strawhits,
		 TrkParticle const&                tpart = _eminus, 
		 TrkFitDirection const&            fdir  = _downstream) : 
      TrkDef(strawcollection,strawhits,tpart,fdir) 
    {
      _shpos  = shposcollection;
      _shfcol = ShFlagCollection; 
    }

    HelixDefHack() {}

    const StrawHitPositionCollection* strawHitPositionCollection() const { return _shpos ; }
    const StrawHitFlagCollection*     strawHitFlagCollection    () const { return _shfcol; }
  };

};

#endif
