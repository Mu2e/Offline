// Filters out the VD101 StepPointMCs ready for resampling
//
// Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

using namespace std;

namespace mu2e{
  class STMResamplingFilter : public art::EDFilter
  {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    struct Config
    {
      fhicl::Atom<std::string> stepsTag{Name("VD101StepPointMCsTag"), Comment("Input tag of StepPointMCs associated with VD101")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("Verbosity of output")};
    };

    using Parameters=art::EDFilter::Table<Config>;

    explicit STMResamplingFilter(const Parameters& pset);
    virtual bool filter(art::Event& event) override;

  private:
    art::InputTag _stepsTag;
  };
  // ===================================================
  STMResamplingFilter::STMResamplingFilter(const Parameters& config) :
    art::EDFilter{config},
    _stepsTag(config().stepsTag())
    {};
  // ===================================================
  bool STMResamplingFilter::filter(art::Event& event)
  {
    // Define a handle to the virtualdetector
    art::Handle<StepPointMCCollection> _inputStepPointMCs;
    event.getByLabel(_stepsTag, _inputStepPointMCs);

    // Check if handle is valid
    if (!(_inputStepPointMCs.isValid()))
    {
      std::cout << _stepsTag << " is an invalid StepPointMC tag, exiting." << std::endl;
      exit(0);
    }

    // Only keep events that have non-zero size
    if(_inputStepPointMCs->size() > 0){return true;}
    else{return false;};
  };
  // ===================================================
}

DEFINE_ART_MODULE(mu2e::STMResamplingFilter)
