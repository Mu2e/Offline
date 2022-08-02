#ifndef Mu2eKinKal_NullHitInfo_hh
#define Mu2eKinKal_NullHitInfo_hh
namespace mu2e {
  // struct describing null ambiguity hit (wire constraint) properties
  struct NullHitInfo {
    NullHitInfo() : toff_(0.0), tvar_(0.0), dvar_(0.0), usetime_(false), useComboDriftTime_(false) {}
    double toff_; // time offset
    double tvar_; // variance on drift time
    double dvar_; // variance on drift distance
    bool usetime_; // use null hit time residual
    bool useComboDriftTime_; // use the combohit drift time
  };
}
#endif
