//

#include "CalPatRec/inc/HelixDefHack.hh"

namespace mu2e {
  

//-----------------------------------------------------------------------------
// default constructor
//-----------------------------------------------------------------------------
  HelixDefHack::HelixDefHack() :  TrkDefHack() {
  }

//-----------------------------------------------------------------------------
// copy contructor
//-----------------------------------------------------------------------------
  HelixDefHack::HelixDefHack(const HelixDefHack& other){
    _shpos  = other.strawHitPositionCollection();
    _shfcol = other.strawHitFlagCollection();
    
    //parameters from TrkDefHack
    _shcol       = other.strawHitCollection();
    _timeCluster = other._timeCluster;
    _h0          = other._h0;
    _tpart       = other._tpart;
    _fdir        = other._fdir;
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
			     const std::vector<StrawHitIndex>& StrawHits       ,
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
