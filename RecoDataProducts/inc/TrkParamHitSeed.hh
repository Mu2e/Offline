//
//  Persistent representation of a KinKal ParameterHit, used in the
//  persistent representation of Kalman Fit
//
#ifndef RecoDataProducts_TrkParamHitSeed_HH
#define RecoDataProducts_TrkParamHitSeed_HH
#include "KinKal/Detector/ParameterHit.hh"
namespace mu2e {
  struct TrkParamHitSeed {
      using PMASK = std::array<bool,KinKal::NParams()>; // parameter mask
    // default constructor: initialization on declaration
    TrkParamHitSeed() {}
    //KinKal constructor
    TrkParamHitSeed(double time, KinKal::Parameters const& params, PMASK const& pmask) :
      time_(time), params_(params), pmask_(pmask) {}

    double time() const { return time_; }
    KinKal::Parameters const& params() const { return params_; }
    PMASK const& pmask() const { return pmask_; }

    //
    //  Payload
    //
    double time_; // time of this constraint: must be supplied on construction and does not change
    KinKal::Parameters params_; // constraint parameters with covariance
    PMASK pmask_; // subset of parmeters to constrain
  };
}
#endif
