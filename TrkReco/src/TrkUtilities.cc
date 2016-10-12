//
//  Collection of tools useful for dealing with various helices
//  Original Author Dave Brown (LBNL) 26 Aug. 2016
//
// Mu2e
#include "TrkReco/inc/TrkUtilities.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "GeneralUtilities/inc/Angles.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
// BTrk
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
//C++
#include <math.h>
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
using namespace std;
namespace mu2e {
  namespace TrkUtilities {

    bool RobustHelix2Traj (RobustHelix const& helix, HepVector& hpvec, float amsign) {
      bool retval(false);
      // compare the input with this configuration's helicity: these must be the same
      // radius = 0 or lambda=0 are degenerate cases that this representation can't handle
      if(helix.radius() > 0.0 && helix.lambda() != 0.0 && hpvec.num_row() == HelixTraj::NHLXPRM) {
	// radius and omega have inverse magnitude, omega is signed by the angular momentum 
	hpvec[HelixTraj::omegaIndex] = amsign/helix.radius();
	// phi0 is the azimuthal angle of the particle velocity vector at the point
	// of closest approach to the origin.  It's sign also depends on the angular
	// momentum.  To translate from the center, we need to reverse coordinates
	hpvec[HelixTraj::phi0Index] = atan2(-amsign*helix.centerx(),amsign*helix.centery());
	// d0 describes the distance to the origin at closest approach.
	// It is signed by the particle angular momentum WRT the origin.
	// The Helix fit radial bias is anti-correlated with d0; correct for it here.
	hpvec[HelixTraj::d0Index] = amsign*(helix.rcent() - helix.radius());
	// the dip angle is measured WRT the perpendicular, signed by the z component of linear momentum
	hpvec[HelixTraj::tanDipIndex] = amsign*helix.lambda()/helix.radius();
	// must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
	// resolve the loop ambiguity such that the POCA is closest to z=0.
	double refphi = helix.fz0()+amsign*M_PI_2;
	double phi = hpvec[HelixTraj::phi0Index];
	double dphi = Angles::deltaPhi(phi,refphi);
	// choose z0 (which loop) so that f=0 is as close to z=0 as possible
	hpvec[HelixTraj::z0Index] = dphi*hpvec[HelixTraj::tanDipIndex]/hpvec[HelixTraj::omegaIndex]; 
	retval = true;
      }
      return retval;
    }

    void RobustHelixFromMom(Hep3Vector const& pos, Hep3Vector const& mom, double charge, double Bz, RobustHelix& helix){
      // speed of light in mm/nsec
      static double clight =299.792;  // This value should come from conditions FIXME!!!	
      // translation factor from MeV/c to curvature radius.
      double momToRad = 1000.0/(charge*Bz*clight);
      // compute some simple useful parameters
      double pt = mom.perp();
      // transverse radius of the helix
      helix._radius = fabs(pt*momToRad);
      //longitudinal wavelength; sign convention goes with angular rotation
      helix._lambda = -mom.z()*momToRad;
      // circle center
      Hep3Vector center = Hep3Vector(pos.x() + mom.y()*momToRad, pos.y() - mom.x()*momToRad, 0.0);
      helix._rcent = center.perp();
      helix._fcent = center.phi();
      // phi at z=0
      double phi = (pos - center).phi() - pos.z()/helix.lambda();
      // reset to be close to 0
      Angles::deltaPhi(phi);
      helix._fz0 = phi;
    }

    void fillSegment(HelixTraj const& htraj, BbrVectorErr const& momerr, KalSegment& kseg) {
      kseg._fmin = htraj.lowRange();
      kseg._fmax = htraj.hiRange();
      kseg._helix = htraj.parameters()->parameter();
      kseg._hcov = htraj.parameters()->covariance();
      kseg._mom = momerr.mag();
      Hep3Vector md = momerr.unit();
      HepVector mdir(3);
      for(size_t icor=0;icor<3;++icor) // CLHEP auto-conversion is broken, FIXME!!
	mdir[icor] = md[icor];

      kseg._momerr = sqrt(momerr.covMatrix().similarity(mdir));
    }
  
    void fillHitSeeds(const KalRep* krep, std::vector<TrkStrawHitSeed>& hitseeds) {
      // extract the TkrStrawHits from the KalRep
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      // loop over the TrkStrawHits and convert them
      for(auto tsh : tshv ) {
	// set the flag according to the status of this hit
	StrawHitFlag hflag;
	if(tsh->isActive())hflag.merge(StrawHitFlag::active);
	if(tsh->poca().status().success())hflag.merge(StrawHitFlag::doca);
	TrkStrawHitSeed seedhit(tsh->index(), tsh->hitT0(), tsh->fltLen(), tsh->hitLen(),
	    tsh->driftRadius(), tsh->poca().doca(), tsh->ambig(),tsh->driftRadiusErr(), hflag);
	hitseeds.push_back(seedhit);
      }
    }
 // compute the overlap between 2 clusters 
    double overlap(TimeCluster const& tc1, TimeCluster const& tc2) {
      double over(0.0);
      double norm = std::min(tc1._strawHitIdxs.size(),tc2._strawHitIdxs.size());
      // count the overlapping hits
      for(auto ih1 = tc1._strawHitIdxs.begin(); ih1 !=tc1._strawHitIdxs.end(); ++ih1){ 
	for(auto ih2 = tc2._strawHitIdxs.begin(); ih2 !=tc2._strawHitIdxs.end(); ++ih2){
	  if(*ih1 == *ih2)
	    over +=1.0;
	    break;
	}
      }
      // add in CaloCluster; count is as much as all the hits
      if(tc1._caloCluster.isNonnull() && tc2._caloCluster.isNonnull()) {
	if(tc1._caloCluster == tc2._caloCluster)
	  over += norm;
	norm *= 2.0;
      }
      return over/norm;
    }
  } // TrkUtilities
}// mu2e
