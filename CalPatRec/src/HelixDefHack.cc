//

#include "CalPatRec/inc/HelixDefHack.hh"

namespace mu2e {
  

//-----------------------------------------------------------------------------
// default constructor
//-----------------------------------------------------------------------------
  HelixDefHack::HelixDefHack() :  TrkDefHack() {
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  HelixDefHack::HelixDefHack(TrkDefHack const& tdef) : 
    TrkDefHack (tdef), 
    _shpos (NULL),
    _shfcol(NULL)
  {
  }


//-----------------------------------------------------------------------------
// real constructor
//-----------------------------------------------------------------------------
  HelixDefHack::HelixDefHack(const StrawHitCollection*         StrawCollection ,
			     const StrawHitPositionCollection* ShposCollection , 
			     const StrawHitFlagCollection*     ShFlagCollection, 
			     const std::vector<hitIndex>&      StrawHits       ,
			     TrkParticle const&                tpart           ,
			     TrkFitDirection const&            fdir            ) : 
    TrkDefHack(StrawCollection,StrawHits,tpart,fdir) 
  {
    _shpos  = ShposCollection;
    _shfcol = ShFlagCollection; 
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  HelixDefHack::~HelixDefHack() {
  }
//-----------------------------------------------------------------------------
// don't care about TrkDefHack - just invalidate pointers to hits and flags
//-----------------------------------------------------------------------------
  void HelixDefHack::init() {
    _shpos  = NULL;
    _shfcol = NULL;
  }

}
