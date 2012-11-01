// Andrei Gaponenko, 2012

#ifndef MCDataProducts_inc_MARSInfo_hh
#define MCDataProducts_inc_MARSInfo_hh

#include <ostream>

namespace mu2e {

  class MARSInfo {
  public:

    double weight() const { return weight_; }

    int protonNumber() const { return protonNumber_; }

    // For MARS inputs we define a subrun as a unit of the "protons on
    // target" simulation.  Data in different subruns of the same run
    // are statistically independent.  Proton numbers are unuque
    // within a subrun.  When reading input file rand1/399_fort.86.gz
    // the "399" is the subrun number.  On the other hand
    // rand1/399_fort.86.gz and rand2/399_fort.86.gz start with the
    // same protons and are not statistically independent.
    int subRunNumber() const { return subRunNumber_; }

    // to distinguish rand1 from rand2 etc.
    int runNumber() const { return runNumber_; }

    MARSInfo(double w, int nj, int sr, int rn) : weight_(w), protonNumber_(nj), subRunNumber_(sr), runNumber_(rn) {}

    // defautl ctr required by ROOT persistency
    MARSInfo() : weight_(), protonNumber_(), subRunNumber_(), runNumber_() {}

  private:
    double weight_;
    int protonNumber_;
    int subRunNumber_;
    int runNumber_;
  };

  inline bool sameProtonAndSimPath(const MARSInfo& a, const MARSInfo& b) {
    return
      (a.protonNumber() == b.protonNumber())&&
      (a.subRunNumber()==b.subRunNumber())&&
      (a.runNumber()==b.runNumber());
  }

  bool partlyCorrelated(const MARSInfo& a, const MARSInfo& b) {
    return
      (a.protonNumber() == b.protonNumber())&&
      (a.subRunNumber()==b.subRunNumber())&&
      (a.runNumber()!=b.runNumber());
  }

  // "less than" comparison that treats correlated protons as equivalent
  struct CmpProtonId {
    bool operator()(const MARSInfo& a, const MARSInfo& b) const {
      return
        (a.subRunNumber() < b.subRunNumber()) || ((a.subRunNumber() == b.subRunNumber()) &&
                                                  (a.protonNumber() < b.protonNumber()));
    }
  };

  // "less than" comparison to distinguish different proton*simulationPath results
  struct CmpProtonIdAndSimPath {
    bool operator()(const MARSInfo& a, const MARSInfo& b) const {
      return
        (a.runNumber() < b.runNumber()) || ((a.runNumber() == b.runNumber()) &&
                                            ((a.subRunNumber() < b.subRunNumber()) || ((a.subRunNumber() == b.subRunNumber()) &&
                                                                                       (a.protonNumber() < b.protonNumber()))));
    }
  };

  std::ostream& operator<<(std::ostream& os, const MARSInfo& mi);
}

#endif/*MCDataProducts_inc_MARSInfo_hh*/
