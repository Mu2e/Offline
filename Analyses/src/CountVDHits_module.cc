// Count number of hits in a VD using StepPointMCs
// TODO - include a parameter to select which VDs to count as a list.
// Pawel Plesniak, 2024

#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <algorithm>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
namespace mu2e {

  //================================================================
  class CountVDHits : public art::EDAnalyzer {
  public:
    explicit CountVDHits(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
    void endJob() override;
  private:
    art::InputTag _stepPointMCsTag;
    std::vector<int> _enabledVDs;
    std::vector<int> _enabledVDStepPointMCsCount;
  };

  //================================================================
  CountVDHits::CountVDHits(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset),
      _stepPointMCsTag(pset.get<std::string>("StepPointMCsTag")),
      _enabledVDs(pset.get<std::vector<int> >("enableVDs", std::vector<int>()))
  {
    // Sort all elements to be in ascending order
    std::sort(_enabledVDs.begin(), _enabledVDs.end());
    // Remove all duplicate elements
    std::vector<int>::iterator _uniqueIndex = std::unique(_enabledVDs.begin(), _enabledVDs.end());
    // Removed undefined memory
    _enabledVDs.resize(std::distance(_enabledVDs.begin(), _uniqueIndex));
    // Insert _enabledVDs.size() zeros into the counter vector
    _enabledVDStepPointMCsCount.insert(_enabledVDStepPointMCsCount.end(), _enabledVDs.size(), 0);
  }

  //================================================================
  void CountVDHits::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    art::Handle<StepPointMCCollection> StepPointMCs;
    event.getByLabel(_stepPointMCsTag, StepPointMCs);
    if (!(StepPointMCs.isValid())){return;}

    // Define the index here so there's less overhead with memory allocation
    int index;

    // Loop over all the StepPointMCs to find which ones cross the VDs of interest. If a StepPointMC does not correspond to a volume ID, it is ignored.
    for (auto &StepPointMC : *StepPointMCs)
    {
      //Find the position off the volume ID
      auto it = find(_enabledVDs.begin(), _enabledVDs.end(), StepPointMC.volumeId());
      // If the volumeID is not in the list, skip to the next event
      if (it == _enabledVDs.end()){continue;}
      // Increment the counter for the corresponding entry
      index = it - _enabledVDs.begin();
      _enabledVDStepPointMCsCount[index]++;
    }
    return;
  }

  //================================================================
  void CountVDHits::endJob()
  {
    std::cout << "Summary " << std::endl;
    for (uint i = 0; i < _enabledVDStepPointMCsCount.size(); i++)
    {
      std::cout << "VD" << _enabledVDs[i] << " - " << _enabledVDStepPointMCsCount[i] << " StepPointMCs" << std::endl;
    }
    return;
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountVDHits)
