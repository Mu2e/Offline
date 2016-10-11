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
#include "RecoDataProducts/inc/StrawHitIndex.hh"

// #include "RecoDataProducts/inc/StrawHitFlag.hh"

#include "CalPatRec/inc/TrkDefHack.hh"

namespace mu2e {
					// add the StrawHitPosition collection to TrkDefHack
  class HelixDefHack : public TrkDefHack {
  private:
    const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _shfcol;

  public:

    HelixDefHack(TrkDefHack const& tdef);

    HelixDefHack(const StrawHitCollection*         strawcollection,
		 const StrawHitPositionCollection* shposcollection, 
		 const StrawHitFlagCollection*     ShFlagCollection, 
		 const std::vector<StrawHitIndex>&      strawhits,
		 TrkParticle const&                tpart = _eminus, 
		 TrkFitDirection const&            fdir  = _downstream);

    HelixDefHack(const HelixDefHack& Copy);

    HelixDefHack();

    ~HelixDefHack();

    void init();

    const StrawHitPositionCollection* strawHitPositionCollection() const { return _shpos ; }
    const StrawHitFlagCollection*     strawHitFlagCollection    () const { return _shfcol; }
    
  };

};

#endif
