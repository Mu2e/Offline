////////////////////////////////////////////////////////////////////////
// Class:       SeparateCaloShowerSteps
// Plugin Type: producer (art v2_06_02)
// File:        SeparateCaloShowerSteps_module.cc
//
// Takes a CaloShowerStepCollections from the output of digi compression
// and separates them out into calorimeter and calorimeterRO instances
//
// A. Edmonds, October 2018
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include <memory>

#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

namespace mu2e {
  class SeparateCaloShowerSteps;
}


class mu2e::SeparateCaloShowerSteps : public art::EDProducer {
public:
  explicit SeparateCaloShowerSteps(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SeparateCaloShowerSteps(SeparateCaloShowerSteps const &) = delete;
  SeparateCaloShowerSteps(SeparateCaloShowerSteps &&) = delete;
  SeparateCaloShowerSteps & operator = (SeparateCaloShowerSteps const &) = delete;
  SeparateCaloShowerSteps & operator = (SeparateCaloShowerSteps &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions

private:

  // art tags for the input collections
  art::InputTag _inputTag;
  std::string _crystalOutputInstance;
  std::string _sipmOutputInstance;


  // unique_ptrs to the new output collections
  std::unique_ptr<CaloShowerStepCollection> _newCrystalSteps;
  std::unique_ptr<CaloShowerStepCollection> _newSiPMSteps;
};


mu2e::SeparateCaloShowerSteps::SeparateCaloShowerSteps(fhicl::ParameterSet const & pset)
  : art::EDProducer{pset},
    _inputTag(pset.get<art::InputTag>("inputTag")),
    _crystalOutputInstance(pset.get<std::string>("crystalOutputInstance")),
    _sipmOutputInstance(pset.get<std::string>("sipmOutputInstance"))
{
  // Call appropriate produces<>() functions here.
  produces<CaloShowerStepCollection>(_crystalOutputInstance);
  produces<CaloShowerStepCollection>(_sipmOutputInstance);
}

void mu2e::SeparateCaloShowerSteps::produce(art::Event & event)
{
  // Implementation of required member function here.
  _newCrystalSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);  
  _newSiPMSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);  

  mu2e::GeomHandle<mu2e::Calorimeter> ch;
  const auto& _calorimeter = ch.get();
  int n_crystals = _calorimeter->nCrystal();

  const auto& caloShowerStepsHandle = event.getValidHandle<CaloShowerStepCollection>(_inputTag);
  const CaloShowerStepCollection& caloShowerSteps = *caloShowerStepsHandle;

  for (const auto& i_caloShowerStep : caloShowerSteps) {
    // if the volume ID is greater than the number of crystals, then we know it is a SiPM step
    CaloShowerStep new_step = i_caloShowerStep;
    int volumeId = i_caloShowerStep.volumeId();
    if (volumeId >= n_crystals) {
      _newSiPMSteps->push_back(new_step);
      continue;
    }

    // sometimes SiPM IDs can be the same as crystal IDs so check that if this step is within the crystal
    CLHEP::Hep3Vector crystalSize = _calorimeter->crystal(volumeId).size(); // size is in full-lengths

    CLHEP::Hep3Vector step_position = i_caloShowerStep.position();
    double delta = 1e-5;
    if (step_position.x() >= -crystalSize.x()/2.-delta && step_position.x() <= crystalSize.x()/2.+delta
	&& step_position.y() >= -crystalSize.y()/2.-delta && step_position.y() <= crystalSize.y()/2.+delta
	&& step_position.z() >= -delta && step_position.z() <= crystalSize.z()+delta) {
      _newCrystalSteps->push_back(new_step);
    }
    else {
      _newSiPMSteps->push_back(new_step);
    }
  }
  
  event.put(std::move(_newCrystalSteps), _crystalOutputInstance);
  event.put(std::move(_newSiPMSteps), _sipmOutputInstance);
}



DEFINE_ART_MODULE(mu2e::SeparateCaloShowerSteps)
