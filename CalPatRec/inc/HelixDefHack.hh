//
// Object to perform helix fit to straw hits
//
// $Id: HelixDefHack.hh,v 1.1 2014/04/04 21:23:34 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/04 21:23:34 $
//
#ifndef HelixDefHack_HH
#define HelixDefHack_HH

// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

// #include "RecoDataProducts/inc/StrawHitFlag.hh"

#include "KalmanTests/inc/TrkDef.hh"

namespace mu2e {
					// add the StrawHitPosition collection to TrkDef
  class HelixDefHack : public TrkDef {
  public:
    HelixDefHack(TrkDef const& tdef) : TrkDef(tdef), _shpos(0) {}

    HelixDefHack(const StrawHitCollection*         strawcollection,
		 const StrawHitPositionCollection* shposcollection, 
		 const std::vector<hitIndex>&      strawhits,
		 TrkParticle const&                tpart = _eminus, 
		 TrkFitDirection const&            fdir  = _downstream) : 
      TrkDef(strawcollection,strawhits,tpart,fdir), 
      _shpos(shposcollection) {}

    HelixDefHack() {}

    const StrawHitPositionCollection* strawHitPositionCollection() const { return _shpos; }
  private:
    const StrawHitPositionCollection* _shpos;
  };

};

#endif
