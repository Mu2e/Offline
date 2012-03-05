#ifndef FilterEfficiency_hh
#define FilterEfficiency_hh

#include <ostream>

namespace mu2e {

  class FilterEfficiency {

    unsigned accepted_;
    unsigned total_;

  public:

    FilterEfficiency() : accepted_(0), total_(0) {}
    FilterEfficiency(unsigned accepted, unsigned total) : accepted_(accepted), total_(total) {}

    unsigned accepted() const { return accepted_; }
    unsigned total() const { return total_; }

    void fill(bool passed) { ++total_; if(passed) ++accepted_; }
  };

  std::ostream& operator<<(std::ostream& os, const FilterEfficiency& eff);

} // namespace mu2e

#endif /* FilterEfficiency_hh */
