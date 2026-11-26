#ifndef Mu2eKinKal_KKFitUtilities_hh
#define Mu2eKinKal_KKFitUtilities_hh
//
//  untemplated utiltity classes and functions
//
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  class ComboHit;
  class Straw;
  class StrawResponse;
  namespace Mu2eKinKal{
    enum Dimension { dresid=0, tresid=1, lresid=2};  // residual dimensions
    // function to turn a StrawHit into a SensorLine object
    std::shared_ptr<KinKal::SensorLine> hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse);
    // function to turn a Straw into a SensorLine object, used to compute material interactions
    std::shared_ptr<KinKal::SensorLine> strawLine(Straw const& straw,double time);
    // test whether a point is inside the detector
    bool inDetector(KinKal::VEC3 const& point);
    // returns phi angle from 0 to pi/2, where 0 is in B field direction
    double LorentzAngle(KinKal::ClosestApproachData const& ptca, KinKal::VEC3 const& bdir);

    bool insideStraw(KinKal::ClosestApproachData const& tpdata,Straw const& straw,double tolerance=0.0);
    // return the time range bounding a set of hits
    KinKal::TimeRange timeBounds(ComboHitCollection const& chits);
    // return the middle of a set of hits
    double zMid(ComboHitCollection const& chits);

  }
}
#endif
