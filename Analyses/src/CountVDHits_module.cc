// Count number of hits in a VD using StepPointMCs
// TODO - include a parameter to select which VDs to count as a list.
// Pawel Plesniak, 2024

#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <algorithm>
#include <limits>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "canvas/Utilities/InputTag.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
namespace mu2e {

  //================================================================
  class CountVDHits : public art::EDAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> stepPointMCsTag{Name("StepPointMCsTag"), Comment("Name of StepPointMCs to search, likely to be g4run:virtualdetector")};
      fhicl::Sequence<int> enabledVDs{Name("enableVDs"), Comment("Vector of which virtual detectors to counnt hits for by number.")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("If verbose, this will print the counting defined in endJob")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CountVDHits(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
    Config _conf;
    art::ProductToken<StepPointMCCollection> _stepPointMCsToken;
    std::vector<int> _enabledVDs;
    std::vector<int> _enabledVDStepPointMCsCount;
    bool _verbose = true;
  };

  //================================================================
  CountVDHits::CountVDHits(const Parameters& conf)
    : art::EDAnalyzer(conf),
      _conf(conf()),
      _stepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
      _enabledVDs(conf().enabledVDs())
  {
    // Sort all elements to be in ascending order
    std::sort(_enabledVDs.begin(), _enabledVDs.end());
    // Remove all duplicate elements
    std::vector<int>::iterator _uniqueIndex = std::unique(_enabledVDs.begin(), _enabledVDs.end());
    // Removed undefined memory
    _enabledVDs.resize(std::distance(_enabledVDs.begin(), _uniqueIndex));
    // Insert _enabledVDs.size() zeros into the counter vector
    _enabledVDStepPointMCsCount.insert(_enabledVDStepPointMCsCount.end(), _enabledVDs.size(), 0);

    // If verbose has been provided, override the local variable
    if (conf().verbose.hasValue()){
      _verbose = *std::move(conf().verbose());
    };
  }

  //================================================================
  void CountVDHits::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCs = event.getProduct(_stepPointMCsToken);

    // Loop over all the StepPointMCs to find which ones cross the VDs of interest. If a StepPointMC does not correspond to a volume ID, it is ignored.
    for (auto const &StepPointMC : StepPointMCs)
    {
      //Find the index of the volume ID
      auto it = find(_enabledVDs.begin(), _enabledVDs.end(), StepPointMC.volumeId());
      // If the volumeID is not in the list, skip to the next event
      if (it == _enabledVDs.end()){continue;}
      // Increment the counter for the corresponding entry
      const int index = it - _enabledVDs.begin();
      _enabledVDStepPointMCsCount[index]++;
    }
    return;
  }

  //================================================================
  void CountVDHits::endJob()
  {
    if (_verbose){
      std::cout << "################### CountVDHits Summary ###################" << std::endl;
      for (uint i = 0; i < _enabledVDStepPointMCsCount.size(); i++)
      {
        std::cout << "VD" << _enabledVDs[i] << " - " << _enabledVDStepPointMCsCount[i] << " StepPointMCs" << std::endl;
      };
      std::cout << "###########################################################" << std::endl;
    };
    return;
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountVDHits)
