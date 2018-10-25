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
    
      int           nHitsLoop(0), nHitsLoopChecked(0);

      //initialize
      _meanHitRadialDist = 0.;
      _nLoops            = 0;
      _nHitsLoopFailed   = 0;

      const ComboHit*     hit(0);

      float         z_first_hit(0), z_last_hit(0), counter(0);
      bool          isFirst(true);
      float         half_pitch  =  M_PI*fabs(_hel->_helix._lambda);
      float         dz_min_toll = 600.;
      unsigned      nhits = _hel->_hhits.size();

      static const XYZVec zdir(0.0,0.0,1.0);
      float          rpullScaleF(1);
      float          cradres(20.);
      float          cperpres(20.);

      //initialize to 0 the chi2 values
      _chi2dXY   = 0.;
      _chi2dZPhi = 0.;

      const mu2e::RobustHelix  *robustHel = &Helix->helix();
      XYZVec                    wdir(0,0,0);
      static XYZVec             zaxis(0.0,0.0,1.0); // unit in z direction

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

	wdir      = hit->wdir();

	//calculate the residuals on the XY and ZPhi planes
	XYZVec cvec  = PerpVector(hit->pos() - robustHel->center(),Geom::ZDir()); // direction from the circle center to the hit
	XYZVec cdir  = cvec.Unit(); // direction from the circle center to the hit
	float  rwdot = wdir.Dot(cdir); // compare directions of radius and wire
	float  dr    = sqrtf(cvec.mag2()) - robustHel->radius();

	float rwdot2 = rwdot*rwdot;
	// compute radial difference and pull
	float werr   = hit->posRes(mu2e::ComboHit::wire);
	float terr   = hit->posRes(mu2e::ComboHit::trans);
	// the resolution is dominated the resolution along the wire
	//      float rres   = std::max(sqrtf(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2)),minrerr);
	float rres   = sqrtf(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2));
	float rpull  = fabs(dr/rres)*rpullScaleF;

	_chi2dXY  += rpull*rpull;

	//RobustHelix: Z-Phi 
	XYZVec wtdir = zaxis.Cross(wdir);   // transverse direction to the wire
	// XYZVec cvec = PerpVector(hit->pos() - helix.center(),Geom::ZDir()); // direction from the circle center to the hit
	// XYZVec cdir = cvec.Unit();          // direction from the circle center to the hit
	XYZVec cperp = zaxis.Cross(cdir);   // direction perp to the radius

	XYZVec hpos = hit->pos(); // this sets the z position to the hit z
	robustHel->position(hpos);                // this computes the helix expectation at that z
	XYZVec dh = hit->pos() - hpos;   // this is the vector between them
	float dtrans = fabs(dh.Dot(wtdir)); // transverse projection
	float dwire = fabs(dh.Dot(wdir));   // projection along wire direction

	// compute the total resolution including hit and helix parameters first along the wire
	float wres2 = std::pow(hit->posRes(mu2e::ComboHit::wire),(int)2) +
	  std::pow(cradres*cdir.Dot(wdir),(int)2) +
	  std::pow(cperpres*cperp.Dot(wdir),(int)2);
	// transverse to the wires
	float wtres2 = std::pow(hit->posRes(mu2e::ComboHit::trans),(int)2) +
	  std::pow(cradres*cdir.Dot(wtdir),(int)2) +
	  std::pow(cperpres*cperp.Dot(wtdir),(int)2);

	_chi2dZPhi += dwire*dwire/wres2 + dtrans*dtrans/wtres2;
	
	
      }//end loop over the hits

      _nStrawHits = nstrawhits;

      if (counter>0){
	_chi2dXY   = _chi2dXY  /counter;
	_chi2dZPhi = _chi2dZPhi/counter;
      }


      if (counter > 0) _meanHitRadialDist /= counter;
      if (nHitsLoop >= _nMinHitsLoop) {
	++_nLoops;
	nHitsLoopChecked +=  nHitsLoop;
      }
      
      _nHitsLoopFailed   = (int)nhits - nHitsLoopChecked;


      //here we evaluate once the impact parameter
      _d0 = robustHel->rcent  () - robustHel->radius ();
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    int    nLoops           () const { return _nLoops;            }
    
    int    nHitsLoopFailed  () const { return _nHitsLoopFailed;   }

    int    nMinHitsLoop     () const { return _nMinHitsLoop;      }
    
    float  meanHitRadialDist() const { return _meanHitRadialDist; }
    
    float  chi2dXY          () const { return _chi2dXY;           }   
    float  chi2dZPhi        () const { return _chi2dZPhi;         }   
    
    float  d0               () const { return _d0;                }

    float  nstrawhits       () const { return _nStrawHits;        }


  private:

    // PDG particle id.
    const HelixSeed* _hel;
    int        _nMinHitsLoop;
    
    int        _nLoops;
    int        _nStrawHits;
    int        _nHitsLoopFailed;
    float      _meanHitRadialDist;
    float      _chi2dXY;
    float      _chi2dZPhi;
    float      _d0;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_HelixTool_hh */
