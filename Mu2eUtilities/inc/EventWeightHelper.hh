// Compute the final event weight by multiplying several EventWeight
// inputs, and produce a diagnostic histogram.
//
// The object is configured at construction time.
// The update() method must be called once at the beginning of
// a new event and before using the weight() call.
//
// Andrei Gaponenko, 2016

#ifndef Mu2eUtilities_inc_EventWeightHelper_hh
#define Mu2eUtilities_inc_EventWeightHelper_hh

#include <string>                       // for string
#include <algorithm>                    // for max
#include <vector>                       // for vector

#include "canvas/Utilities/InputTag.h"  // for InputTag

namespace fhicl {
class ParameterSet;
}  // namespace fhicl
namespace art {
class Event;
}  // namespace art
class TH1;
namespace art {
class TFileDirectory;
}  // namespace art

namespace mu2e {

  class EventWeightHelper {
  public:
    // must be called once per event, before any calls to the weight() method
    void update(const art::Event& evt);

    double weight() const { return weight_; }

    // histograms will be placed in topdir/subdir
    EventWeightHelper(const fhicl::ParameterSet& pset,
                      art::TFileDirectory topdir,
                      const std::string& subdir);

  private:
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    double weight_;

    TH1 *hfinal_;
  };
}

#endif/*Mu2eUtilities_inc_EventWeightHelper_hh*/
