//
//  Collection of tools useful for dealing with various helices
//  Original Author Dave Brown (LBNL) 26 Aug. 2016
//
// Mu2e
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkStraw.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
// KinKal
#include "KinKal/Trajectory/CentralHelix.hh"
// CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
//C++
#include <cmath>
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e {
  namespace TrkUtilities {
    void RobustHelixFromMom(Hep3Vector const& pos, Hep3Vector const& mom, double charge, double Bz, RobustHelix& helix){
      double momToRad = 1000.0/(charge*Bz*CLHEP::c_light);
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
   // compute the overlap between 2 clusters
    double overlap(SHIV const& shiv1, SHIV const& shiv2) {
      double over(0.0);
      double norm = std::min(shiv1.size(),shiv2.size());
      // count the overlapping hits
      for(auto h1 : shiv1){
        for(auto h2 : shiv2){
          if(h1 == h2){
            over +=1.0;
            break;
          }
        }
      }
      return over/norm;
    }

    double overlap(TimeCluster const& tc1, TimeCluster const& tc2) {
      double hover = overlap(tc1._strawHitIdxs,tc2._strawHitIdxs);
      double norm = std::min(tc1.hits().size(),tc2.hits().size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(tc1._caloCluster.isNonnull() && tc2._caloCluster.isNonnull()) {
        if(tc1._caloCluster == tc2._caloCluster)
          over += norm;
        norm *= 2;
      }
      return over/norm;
    }
    double energy(double mass, double momentum) {  return sqrt(momentum*momentum + mass*mass); }
    double beta(double mass, double momentum) { return fabs(momentum)/energy(mass,momentum); }
    double betagamma(double mass, double momentum) { return fabs(momentum)/mass; }
    double gamma(double mass, double momentum) { return energy(mass,momentum)/mass; }
  } // TrkUtilities
}// mu2e

