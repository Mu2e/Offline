#ifndef Mu2eUtilities_HelixTool_hh
#define Mu2eUtilities_HelixTool_hh
//
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Pezzullo
//

// CLHEP includes
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

namespace mu2e {

  class HelixTool{

  public:
    HelixTool(HelixSeed *Helix,
	      int        MinHitsSingleLoop=3): 
      _hel(Helix), _nMinHitsLoop(MinHitsSingleLoop){
    
      int           nHitsLoop(0), nHitsLoopChecked(0);

      //initialize
      _meanHitRadialDist = 0.;
      _nLoops            = 0;
      _nHitsLoopFailed   = 0;

      ComboHit*     hit(0);

      float         z_first_hit(0), z_last_hit(0), counter(0);
      bool          isFirst(true);
      float         half_pitch  =  M_PI*fabs(_hel->_helix._lambda);
      float         dz_min_toll = 600.;
      unsigned      nhits = _hel->_hhits.size();

      for (unsigned f=0; f<nhits; ++f){
	hit = &_hel->_hhits[f];
	if (hit->_flag.hasAnyProperty(StrawHitFlag::outlier))     continue;
      
	_meanHitRadialDist += sqrtf(hit->pos().x()*hit->pos().x() + hit->pos().y()*hit->pos().y());
	++counter;
	float z = hit->pos().z();
	if (isFirst){
	  z_first_hit = z;
	  z_last_hit  = z;
	  nHitsLoop   = 1;
	  isFirst     = false;
	}else {
	  float    dz_last_hit  = z - z_last_hit;
	  float    dz_first_hit = z - z_first_hit;

	  if ( ( dz_first_hit < half_pitch) && ( dz_last_hit < dz_min_toll)){
	    ++nHitsLoop;
	    z_last_hit        = z;
	  } else {
	    if (nHitsLoop >= _nMinHitsLoop) {
	      ++_nLoops;
	      nHitsLoopChecked +=  nHitsLoop;
	    }
	    nHitsLoop = 0;

	    if ( (z - z_last_hit) >= half_pitch){
	      //reset the hits-per-loop counter
	      nHitsLoop = 1;
	      
	      //re-set the values of the first and last z of the hits within the loop
	      z_first_hit = z;
	      z_last_hit  = z;
	    }
	  }
	}
	
      }//end loop over the hits

      if (counter > 0) _meanHitRadialDist /= counter;
      if (nHitsLoop >= _nMinHitsLoop) {
	++_nLoops;
	nHitsLoopChecked +=  nHitsLoop;
      }
      
      _nHitsLoopFailed   = (int)nhits - nHitsLoopChecked;
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    int    nLoops           () const { return _nLoops;            }
    
    int    nHitsLoopFailed  () const { return _nHitsLoopFailed;   }

    int    nMinHitsLoop     () const { return _nMinHitsLoop;      }
    
    float  meanHitRadialDist() const { return _meanHitRadialDist; }
    
  private:

    // PDG particle id.
    HelixSeed* _hel;
    int        _nMinHitsLoop;
    
    int        _nLoops;
    int        _nHitsLoopFailed;
    float      _meanHitRadialDist;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_HelixTool_hh */
