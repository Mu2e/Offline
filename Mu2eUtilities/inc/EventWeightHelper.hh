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

#include <string>

#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileDirectory.h"

namespace fhicl { class ParameterSet; }
namespace art { class Event; }
class TH1;

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
