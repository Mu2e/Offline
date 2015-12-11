//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.hh,v 1.8 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
#ifndef HelixDef_HH
#define HelixDef_HH
// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// base class 
#include "TrkReco/inc/TrkDef.hh"
namespace mu2e 
{
// add the StrawHitPosition collection to TrkDef
  class HelixDef : public TrkDef {
    public:
      HelixDef(TrkDef const& tdef) : TrkDef(tdef), _shpos(0) {}
      HelixDef(const StrawHitCollection* strawcollection,const StrawHitPositionCollection* shposcollection, const std::vector<hitIndex>& strawhits,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream) : TrkDef(strawcollection,strawhits,tpart,fdir), _shpos(shposcollection) {}
      HelixDef() {}
      const StrawHitPositionCollection* strawHitPositionCollection() const { return _shpos; }
    private:
      const StrawHitPositionCollection* _shpos;
  };
}

#endif
