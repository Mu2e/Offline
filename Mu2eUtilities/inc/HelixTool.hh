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

#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;

namespace mu2e {

  class HelixTool{

  public:
    HelixTool(const HelixSeed *Helix,
	      int        MinHitsSingleLoop=3): 
      _hel(Helix), _nMinHitsLoop(MinHitsSingleLoop){
    
      //initialize
      _meanHitRadialDist = 0.;
      _nLoops            = 0;
      const ComboHit*     hit(0);

      float         z_first_hit(0), z_last_hit(0), counter(0);
      bool          isFirst(true);
      unsigned      nhits = _hel->_hhits.size();

      const mu2e::RobustHelix  *robustHel = &Helix->helix();

      int nstrawhits = 0;

      for (unsigned f=0; f<nhits; ++f){
	hit = &_hel->_hhits[f];
	if (hit->_flag.hasAnyProperty(StrawHitFlag::outlier))     continue;
	nstrawhits += hit -> nStrawHits();

	_meanHitRadialDist += sqrtf(hit->pos().x()*hit->pos().x() + hit->pos().y()*hit->pos().y());
	++counter;
	float z = hit->pos().z();
	if (isFirst){
	  z_first_hit = z;
	  z_last_hit  = z;
	  isFirst     = false;
	}else {
	  z_last_hit  = z;
	}
      }//end loop over the hits

      _nStrawHits = nstrawhits;

      if (counter > 0) _meanHitRadialDist /= counter;

      _nLoops = (z_last_hit - z_first_hit)/(fabs(robustHel->lambda())*2.*M_PI);

      //here we evaluate once the impact parameter
      _d0 = robustHel->rcent  () - robustHel->radius ();


      // we now estiamte the ratio of the number of measured hits to the number of the expected ones
      // we make a few assumptions and appriximations:
      //   - 36 tacking faces
      //   - an average spacing of 178 mm
      //   - we approximate the tracker as a cylinder: r_min = 380 mm, r_max = 860 mm
      int      nFaces(36); 
      float    expected_hits(0);
      float    z_step(178.); 
      float    minR_tracker(380.);
      float    maxR_tracker(860.);

      for (int i=0; i<nFaces;i++){
	float   z = z_first_hit + (double)i*z_step;
	if (z < z_first_hit )  continue;
	if (z > z_last_hit  )  continue;

	XYZVec  pos;
	pos.SetZ(z);
	robustHel->position(pos);
      
	//now check that we are in the active area of tracker
	double   hitR = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
	if (hitR < minR_tracker)  continue;
	if (hitR > maxR_tracker)  continue;

	++expected_hits;
      }

      //each face has two planes. We are assuming 2 ComboHits per face
      _hitRatio = nhits / (2.*expected_hits);
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    float  nLoops           () const { return _nLoops;            }
    int    nMinHitsLoop     () const { return _nMinHitsLoop;      }
    float  meanHitRadialDist() const { return _meanHitRadialDist; }
    float  d0               () const { return _d0;                }
    float  nstrawhits       () const { return _nStrawHits;        }

    //function that evaluates the ratio between the measured ComboHits and the
    //expected intesections of the helix with the tracker planes. This function
    //models the tracker as a perfect cylinder
    float  hitRatio         () const { return _hitRatio;          }


  private:
    const HelixSeed* _hel;
    int        _nMinHitsLoop;

    float      _nLoops;
    int        _nStrawHits;
    float      _meanHitRadialDist;
    float      _d0;
    float      _hitRatio;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_HelixTool_hh */
