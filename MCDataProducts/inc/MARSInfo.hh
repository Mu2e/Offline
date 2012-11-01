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

    MARSInfo(double w, int nj, int sr) : weight_(w), protonNumber_(nj), subRunNumber_(sr) {}

    // defautl ctr required by ROOT persistency
    MARSInfo() : weight_(), protonNumber_(), subRunNumber_() {}

  private:
    double weight_;
    int protonNumber_;
    int subRunNumber_;
  };

  std::ostream& operator<<(std::ostream& os, const MARSInfo& mi);
}

#endif/*MCDataProducts_inc_MARSInfo_hh*/
