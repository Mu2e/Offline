//
//  Collection of tools useful for dealing with various tracking functions
//  not specific to a module or task.
//  Original Author Dave Brown (LBNL) 26 Aug. 2016
//
#ifndef TrkReco_TrkUtilities_HH
#define TrkReco_TrkUtilities_HH
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include <vector>

class HelixTraj;
class BbrVectorErr;
class KalRep;
class TrkDifPieceTraj;
namespace mu2e {
  class RobustHelix;
  class KalSegment;
  class KalSeed;
  class HelixSeed;
  class TrkStrawHitSeed;
  class TrkStraw;
  class TimeCluster;
  class TrkCaloHit;
  class TrkCaloHitSeed;
  class ComboHitCollection;
  typedef std::vector<StrawHitIndex> SHIV;
  namespace TrkUtilities {
    // convert the robust helix format into the BaBar format HelixTraj.  This requires
    // the sign of the angular momentum about the z axis, as the BaBar rep os semi-kinematic
    bool RobustHelix2Traj (RobustHelix const& helix, CLHEP::HepVector& hpvec, float amsign);
    // create a robust helix from raw particle informaiton.  This is useful for MC comparison
    void RobustHelixFromMom(CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& mom, double charge, double Bz, RobustHelix& helix);
    // create a KalSegment (helix segment) from a HelixTraj
    //void fillSegment(HelixTraj const& htraj, double locflt, double globflt, TrkT0 t0, double mass, int charge, BField const& bfield,  KalSegment& kseg);
    // create HitSeeds from the TrkStrawHits in a KalRep
    void fillStrawHitSeeds(const KalRep* krep, ComboHitCollection const& chits, std::vector<TrkStrawHitSeed>& hitseeds);
    void fillCaloHitSeed(const TrkCaloHit* chit, CLHEP::Hep3Vector const& tmom, TrkCaloHitSeed& caloseed);
    void fillStraws(const KalRep* krep, std::vector<TrkStraw>& straws);
    // compute overlap between 2 time clusters
    double overlap(TimeCluster const& tc1, TimeCluster const& tc2);
    double overlap(KalSeed const& ks1, KalSeed const& ks2);
    double overlap(KalSeed const& ks, HelixSeed const& hs);
    double overlap(HelixSeed const& hs,TimeCluster const& tc);
    double overlap(SHIV const& shiv1, SHIV const& shiv2);
    // compute the flightlength for a given z position
    //    double zFlight(TrkDifPieceTraj const& ptraj, double pz);
    double chisqConsistency(const KalRep* krep);
    unsigned countBends(const KalRep* krep);
    const TrkCaloHit* findTrkCaloHit(const KalRep* krep);
    // simple kinematic utilities
    double energy(double mass, double momentum);
    double beta(double mass, double momentum);
    double betagamma(double mass, double momentum);
    double gamma(double mass, double momentum);
  }
}
#endif
