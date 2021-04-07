//Author: S Middleton
//Purpose: Reco product for the cosmic kalman
//TODO - do we need this? Maybe we can re-use KalSeed 
#ifndef RecoDataProducts_CosmicKalSeed_HH
#define RecoDataProducts_CosmicKalSeed_HH
// mu2e
#include "TrkReco/inc/TrkPrintUtils.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "RecoDataProducts/inc/TrkStraw.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/KalSegment.hh"
#include "canvas/Persistency/Common/Ptr.h"
// BTrk
#include "DataProducts/inc/PDGCode.hh"
// Root
#include <Rtypes.h>
// C++
#include <vector>
namespace mu2e {

  struct CosmicKalSeed {
	  CosmicKalSeed() :  _chisq(-1.0), _fitcon(-1.0), _flt0(0)  {}
    CosmicKalSeed(PDGCode::type tpart,TrkFitDirection fdir, TrkFitFlag const& status, double flt0=0.0 ) :
      _tpart(tpart), _fdir(fdir), _status(status),
      _chisq(-1.0), _fitcon(-1.0), _flt0(static_cast<Float_t>(flt0)){}


	  PDGCode::type const& particle() const { return _tpart; }
	  TrkFitDirection const& fitDirection() const { return _fdir; }
	  std::vector<TrkStrawHitSeed> const& hits() const { return _hits;}
	  ComboHitCollection strawcombohits() const { return _StrawLevelComboHits;}
	  std::vector<TrkStraw> const& straws() const { return _straws;}
	  std::vector<KalSegment> const& segments() const { return _segments; }
	  //std::vector<KalSegment>::const_iterator nearestSegment(float fltlen)  const;
	  //std::vector<KalSegment>::const_iterator nearestSegment(const XYZVec& pos)  const; // find nearest segment to a GLOBAL position
	  TrkFitFlag const& status() const { return _status; }
	  Float_t flt0() const { return _flt0; }
	  HitT0 t0() const;
	  Float_t chisquared() const { return _chisq; }
	  Float_t fitConsistency() const { return _fitcon; }

	  art::Ptr<CosmicTrackSeed> const& cosmicseed() const { return _cosmicseed; }

	  // global information about the track
	  PDGCode::type			 _tpart; // particle assumed for this fit
	  TrkFitDirection	 _fdir; // direction in which this particle was fit TODO
	  TrkFitFlag			 _status; // status of this fit
	  HitT0			    	 _t0; // track t0; Time particle crosses z=0
	  Float_t			    	_chisq; // fit chisquared value
	  Float_t			    _fitcon; // fit consistency
    Float_t			    	_flt0; // flight distance where the track crosses the tracker midplane (z=0)
    
	  std::vector<KalSegment>	    	_segments; // segments of the Kalman filter fit result
	  std::vector<TrkStrawHitSeed>    _hits; // hit seeds for all the hits used in this fit
	  std::vector<TrkStraw>	    	_straws; // straws interesected by this fit

	  art::Ptr<CosmicTrackSeed>       _cosmicseed; // associated Cosmic Seed (for seed fits); can be null
	  
	  ComboHitCollection		_StrawLevelComboHits;

    bool kalmanworked = false;
	};
	typedef std::vector<mu2e::CosmicKalSeed> CosmicKalSeedCollection;
}

#endif /*CosmicKalSeed.hh*/
