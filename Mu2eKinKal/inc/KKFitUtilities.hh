#ifndef Mu2eKinKal_KKFitUtilities_hh
#define Mu2eKinKal_KKFitUtilities_hh
//
//  untemplated utiltity classes and functions
//
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  class ComboHit;
  class Straw;
  class StrawResponse;
  namespace Mu2eKinKal{
    enum Dimension { dresid=0, tresid=1};  // residual dimensions
    // function to turn a StrawHit into a Line object
    KinKal::Line hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse);
    // test whether a point is inside the detector
    bool inDetector(KinKal::VEC3 const& point);
    double LorentzAngle(KinKal::ClosestApproachData const& ptca, KinKal::VEC3 const& bdir);

    // this finds the time at which the traj crosses the given z value
    template <class KTRAJ> double zTime(KTRAJ const& ktraj, double zpos, double thint) {
      double zhint = ktraj.position3(thint).Z();
      double vz = ktraj.velocity(thint).Z();
      double ztime = thint + (zpos-zhint)/vz;
      return ztime;
    }
    // this finds the time at which the traj crosses the given z value
    template <class KTRAJ> double zTime(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, double zpos, double thint) {
      double ztime(thint);
      auto zindex = ptraj.nearestIndex(ztime);
      size_t oldzindex;
      unsigned ntries(0);
      static const unsigned maxntries(100);// this usually converges in 1 iteration, but it can oscillate with bad fits
      do {
        oldzindex = zindex;
        auto const& traj = ptraj.piece(zindex);
        ztime = zTime(traj,zpos,ztime);
        zindex = ptraj.nearestIndex(ztime);
        ++ntries;
      } while (ntries < maxntries && zindex != oldzindex);
//      if(ntries == maxntries) throw cet::exception("RECO")<<"mu2e::KKFitUtilities: zTime failure" << std::endl;
      return ztime;
    }
    bool insideStraw(KinKal::ClosestApproachData const& tpdata,Straw const& straw,double tolerance=0.0);

  }
}
#endif
