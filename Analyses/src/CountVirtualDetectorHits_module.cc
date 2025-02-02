// Count number of hits in a VD using StepPointMCs
// Original author: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <bits/stdc++.h>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"


namespace mu2e {
  class CountVirtualDetectorHits : public art::EDAnalyzer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> stepPointMCsTag{Name("StepPointMCsTag"), Comment("Name of StepPointMCs to search, likely to be g4run:virtualdetector")};
      fhicl::Sequence<int> enabledVDs{Name("virtualDetectorIDs"), Comment("Vector of which virtual detectors to counnt hits for by number.")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CountVirtualDetectorHits(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
    Config _conf;
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    std::vector<int> enabledVDs, counter;
    int index = 0, nVDs = 0;
  };

  CountVirtualDetectorHits::CountVirtualDetectorHits(const Parameters& conf)
    : art::EDAnalyzer(conf),
      StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
      enabledVDs(conf().enabledVDs()) {
        // Create list of unique enabled virtual detectors
        std::set<int> enabledVDsSet(enabledVDs.begin(), enabledVDs.end());
        enabledVDs.clear();
        enabledVDs.insert(enabledVDs.end(), enabledVDsSet.begin(), enabledVDsSet.end());

        // Insert _enabledVDs.size() zeros into the counter vector
        nVDs = enabledVDs.size();
        counter.clear();
        counter.insert(counter.end(), nVDs, 0);
      };

  void CountVirtualDetectorHits::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.size() == 0)
      throw cet::exception("DataError") << "Requested data not found";

    // Loop over all the StepPointMCs to find which ones cross the VDs of interest. If a StepPointMC does not correspond to a volume ID, it is ignored.
    for (auto const &StepPointMC : StepPointMCs) {
      index = std::find(enabledVDs.begin(), enabledVDs.end(), StepPointMC.volumeId()) - enabledVDs.begin();
      // Discard the entry if it is not an entry in the list of enabledVDs
      if (index == nVDs)
        continue;
      // Increment the counter for the corresponding entry
      counter[index]++;
    };
    return;
  };

  void CountVirtualDetectorHits::endJob() {
    mf::LogInfo log("Unfiltered virtual detector hits summary");
    log << "\nNumber of hits in sim by virtual detector ID\n";
    for (int i = 0; i < nVDs; i++)
      log << "VD" << enabledVDs[i] << ": " << counter[i] << "\n";
    log << "\n";
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountVirtualDetectorHits)
