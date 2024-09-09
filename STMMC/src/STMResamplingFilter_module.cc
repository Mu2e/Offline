// Removes empty StepPointMC collections with the associated tag
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
      fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("Prints summary")};
    };

    using Parameters=art::EDFilter::Table<Config>;

    explicit STMResamplingFilter(const Parameters& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  private:
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    bool verbose = false;
    uint keptEvents = 0;
    uint discardedEvents = 0;
  };
  // ===================================================
  STMResamplingFilter::STMResamplingFilter(const Parameters& conf) :
    art::EDFilter{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag()))
    {
      auto _verbose = conf().verbose();
      if(_verbose)verbose = *_verbose;
    };
  // ===================================================
  bool STMResamplingFilter::filter(art::Event& event)
  {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);

    // Only keep events that have non-zero size
    if(StepPointMCs.size() > 0){
      keptEvents++;
      return true;
    }
    else{
      discardedEvents++;
      return false;
    };
  };
  // ===================================================
  void STMResamplingFilter::endJob()
  {
    if (verbose == true)
      {
        std::cout << "Number of kept events: " << keptEvents << std::endl;
        std::cout << "Number of discarded events: " << discardedEvents << std::endl;
      };
  };
}

DEFINE_ART_MODULE(mu2e::STMResamplingFilter)
