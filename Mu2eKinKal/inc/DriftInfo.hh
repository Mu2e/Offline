#ifndef Mu2eKinKal_DriftInfo_hh
#define Mu2eKinKal_DriftInfo_hh
namespace mu2e {
// struct describing local drift info
  struct DriftInfo {
    DriftInfo() : tdrift_(0.0), tdriftvar_(0.0), vdrift_(0.0) {}
    double tdrift_; // drift time
    double tdriftvar_; // variance on drift time
    double vdrift_; // instantanious drift speed
  };
}
#endif
