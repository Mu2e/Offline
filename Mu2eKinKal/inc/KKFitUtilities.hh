#ifndef Mu2eKinKal_KKFitUtilities_hh
#define Mu2eKinKal_KKFitUtilities_hh
//
//  untemplated utiltity classes and functions
//
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
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
  }
}
#endif
