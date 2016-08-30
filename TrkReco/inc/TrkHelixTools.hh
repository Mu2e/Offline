//
//  Collection of tools useful for dealing with various helices
//  Original Author Dave Brown (LBNL) 26 Aug. 2016
//
#ifndef TrkReco_TrkHelixTools_HH
#define TrkReco_TrkHelixTools_HH
#include "CLHEP/Vector/ThreeVector.h"

class HelixTraj;
namespace mu2e {
  class RobustHelix;
  namespace TrkHelixTools {
  // convert the robust helix format into the BaBar format HelixTraj.  This requires
  // the sign of the angular momentum about the z axis, as the BaBar rep os semi-kinematic
    bool RobustHelix2Traj (RobustHelix const& helix, HelixTraj &traj, float amsign);
    void RobustHelixFromMom(CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& mom, double charge, double Bz, RobustHelix& helix);
  }
}
#endif
