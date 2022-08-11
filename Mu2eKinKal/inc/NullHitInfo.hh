#ifndef Mu2eKinKal_NullHitInfo_hh
#define Mu2eKinKal_NullHitInfo_hh
namespace mu2e {
  // struct describing null ambiguity hit (wire constraint) properties
  struct NullHitInfo {
    enum nullTimeMode{none=0,usedoca,usecombo};
   NullHitInfo() : toff_(0.0), tvar_(0.0), dvar_(0.0), tmode_(none) {}
    NullHitInfo(double toff, double tvar, double dvar, nullTimeMode ntm) :
      toff_(toff), tvar_(tvar), dvar_(dvar), tmode_(ntm) {}
    double toff_; // time offset
    double tvar_; // variance on drift time
    double dvar_; // variance on drift distance
    nullTimeMode tmode_; // mode for using time
    bool useTime() const { return tmode_ > none; }
    bool useComboDriftTime() const { return tmode_ == usecombo; }
  };
}
#endif
