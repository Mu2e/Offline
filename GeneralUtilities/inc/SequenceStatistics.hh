#ifndef GeneralUtilities_SequenceStatistics_hh
#define GeneralUtilities_SequenceStatistics_hh

//
// Compute min/max and mean/RMS about a sequence of numbers.
// In the future, this should be extended to include general moments.
//
//
// Contact person Rob Kutschke
//

#include "Offline/GeneralUtilities/inc/MinMax.hh"
#include "Offline/GeneralUtilities/inc/RMS.hh"

#include <iostream>

class SequenceStatistics{

public:

  SequenceStatistics():
    rms_(),
    minmax_()
  {}

  SequenceStatistics(double x):
    rms_(x),
    minmax_(x)
  {}

  // Accept compiler written d'tor, copy c'tor and assignment operator.

  // Forward this function since both parts must give the same answer.
  int n() const { return rms_.n(); }

  // Otherwise just return the underlying objects and don't bother forwarding.
  RMS    const& moments() const { return rms_;}
  MinMax const& limits()  const { return minmax_;}

  void accumulate(double x){
    rms_.accumulate(x);
    minmax_.accumulate(x);
  }

private:
  RMS rms_;
  MinMax minmax_;

};

inline std::ostream& operator<<(std::ostream& ost,
                               const SequenceStatistics& stats ){
  ost << stats.moments() << " "
      << stats.limits();
  return ost;
}

#endif /* GeneralUtilities_SequenceStatistics_hh */
