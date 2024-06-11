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
      fhicl::Atom<art::InputTag> stepPointMCsTag{Name("VD101StepPointMCsTag"), Comment("Input tag of StepPointMCs associated with VD101")};
    };

    using Parameters=art::EDFilter::Table<Config>;

    explicit STMResamplingFilter(const Parameters& pset);
    virtual bool filter(art::Event& event) override;

  private:
    art::ProductToken<StepPointMCCollection> _stepPointMCsToken;
  };
  // ===================================================
  STMResamplingFilter::STMResamplingFilter(const Parameters& conf) :
    art::EDFilter{conf},
    _stepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag()))
    {};
  // ===================================================
  bool STMResamplingFilter::filter(art::Event& event)
  {
    auto const& StepPointMCs = event.getProduct(_stepPointMCsToken);

    // Only keep events that have non-zero size
    if(StepPointMCs.size() > 0){return true;}
    else{return false;};
  };
  // ===================================================
}

DEFINE_ART_MODULE(mu2e::STMResamplingFilter)
