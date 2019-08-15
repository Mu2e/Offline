////////////////////////////////////////////////////////////////////////
// Class:       FixCaloShowerStepPtrs
// Plugin Type: producer (art v2_06_02)
// File:        FixCaloShowerStepPtrs_module.cc
//
// This module fixes an off-by-one bug in earlier versions of CompressDigiMCs
// before commit XXXXXXX
//
// Generated at Wed Apr 12 16:10:46 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
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

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "DataProducts/inc/IndexMap.hh"

namespace mu2e {
  class FixCaloShowerStepPtrs;

  typedef std::map<art::Ptr<mu2e::CaloShowerStep>, art::Ptr<mu2e::CaloShowerStep> > CaloShowerStepRemap;
}


class mu2e::FixCaloShowerStepPtrs : public art::EDProducer {
public:
  explicit FixCaloShowerStepPtrs(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FixCaloShowerStepPtrs(FixCaloShowerStepPtrs const &) = delete;
  FixCaloShowerStepPtrs(FixCaloShowerStepPtrs &&) = delete;
  FixCaloShowerStepPtrs & operator = (FixCaloShowerStepPtrs const &) = delete;
  FixCaloShowerStepPtrs & operator = (FixCaloShowerStepPtrs &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  art::Ptr<mu2e::CaloShowerStep> copyCaloShowerStep(const mu2e::CaloShowerStep& old_calo_shower_step);
  void copyCaloShowerSim(const mu2e::CaloShowerSim& old_calo_shower_sim, const CaloShowerStepRemap& remap);
  void copyCaloShowerStepRO(const mu2e::CaloShowerStepRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap);

private:

  // art tags for the input collections
  std::vector<art::InputTag> _caloShowerStepTags;
  art::InputTag _caloShowerSimTag;
  art::InputTag _caloShowerStepROTag;

  // handles to the old collections
  art::Handle<CaloShowerStepCollection> _caloShowerStepsHandle;
  art::Handle<CaloShowerSimCollection> _caloShowerSimsHandle;
  art::Handle<CaloShowerStepROCollection> _caloShowerStepROsHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CaloShowerSimCollection> _newCaloShowerSims;
  std::unique_ptr<CaloShowerStepROCollection> _newCaloShowerStepROs;

  // for StepPointMCs, SimParticles and GenParticles we also need reference their new locations with art::Ptrs and so need their ProductIDs and Getters
  art::ProductID _newCaloShowerStepsPID;
  const art::EDProductGetter* _newCaloShowerStepGetter;
  std::map<art::ProductID, const art::EDProductGetter*> _oldCaloShowerStepGetter;
};


mu2e::FixCaloShowerStepPtrs::FixCaloShowerStepPtrs(fhicl::ParameterSet const & pset)
  : art::EDProducer{pset},
    _caloShowerStepTags(pset.get<std::vector<art::InputTag> >("caloShowerStepTags")),
    _caloShowerSimTag(pset.get<art::InputTag>("caloShowerSimTag")),
    _caloShowerStepROTag(pset.get<art::InputTag>("caloShowerStepROTag"))
{
  // Call appropriate produces<>() functions here.
  produces<CaloShowerStepCollection>();
  produces<CaloShowerSimCollection>();
  produces<CaloShowerStepROCollection>();
}

void mu2e::FixCaloShowerStepPtrs::produce(art::Event & event)
{

  // Check to see if this collection has the bug
  _newCaloShowerSims = std::unique_ptr<CaloShowerSimCollection>(new CaloShowerSimCollection);
  event.getByLabel(_caloShowerSimTag, _caloShowerSimsHandle);
  const auto& caloShowerSims = *_caloShowerSimsHandle;
  bool bug_found = true;
  for (const auto& i_caloShowerSim : caloShowerSims) {
    for (const auto& i_caloShowerStep : i_caloShowerSim.caloShowerSteps()) {
      if(i_caloShowerStep.key()==0) {
	bug_found = false;
	break;
      }
    }
  }


  // Implementation of required member function here.
  CaloShowerStepRemap caloShowerStepRemap;
  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  _newCaloShowerStepsPID = event.getProductID<CaloShowerStepCollection>();
  _newCaloShowerStepGetter = event.productGetter(_newCaloShowerStepsPID);
  for (std::vector<art::InputTag>::const_iterator i_tag = _caloShowerStepTags.begin(); i_tag != _caloShowerStepTags.end(); ++i_tag) {
    const auto& oldCaloShowerSteps = event.getValidHandle<CaloShowerStepCollection>(*i_tag);
    art::ProductID i_product_id = oldCaloShowerSteps.id();
    _oldCaloShowerStepGetter[i_product_id] = event.productGetter(i_product_id);
      
    for (CaloShowerStepCollection::const_iterator i_caloShowerStep = oldCaloShowerSteps->begin(); i_caloShowerStep != oldCaloShowerSteps->end(); ++i_caloShowerStep) {
      if (bug_found) {
	art::Ptr<mu2e::CaloShowerStep> oldShowerStepPtr(i_product_id,  i_caloShowerStep - oldCaloShowerSteps->begin()+1, _oldCaloShowerStepGetter[i_product_id]);
	art::Ptr<mu2e::CaloShowerStep> newShowerStepPtr = copyCaloShowerStep(*i_caloShowerStep);
	caloShowerStepRemap[oldShowerStepPtr] = newShowerStepPtr;
      }
      else {
	art::Ptr<mu2e::CaloShowerStep> oldShowerStepPtr(i_product_id,  i_caloShowerStep - oldCaloShowerSteps->begin(), _oldCaloShowerStepGetter[i_product_id]);
	art::Ptr<mu2e::CaloShowerStep> newShowerStepPtr = copyCaloShowerStep(*i_caloShowerStep);
	caloShowerStepRemap[oldShowerStepPtr] = newShowerStepPtr;
	//	std::cout << "AE: " << oldShowerStepPtr << " --> " << newShowerStepPtr << std::endl;
      }
    }
  }    

  for (const auto& i_caloShowerSim : caloShowerSims) {
    copyCaloShowerSim(i_caloShowerSim, caloShowerStepRemap);
  }

  _newCaloShowerStepROs = std::unique_ptr<CaloShowerStepROCollection>(new CaloShowerStepROCollection);
  event.getByLabel(_caloShowerStepROTag, _caloShowerStepROsHandle);
  const auto& caloShowerStepROs = *_caloShowerStepROsHandle;
  for (const auto& i_caloShowerStepRO : caloShowerStepROs) {
    copyCaloShowerStepRO(i_caloShowerStepRO, caloShowerStepRemap);
  }
    
  event.put(std::move(_newCaloShowerSteps));
  event.put(std::move(_newCaloShowerSims));
  event.put(std::move(_newCaloShowerStepROs));
}

art::Ptr<mu2e::CaloShowerStep> mu2e::FixCaloShowerStepPtrs::copyCaloShowerStep(const mu2e::CaloShowerStep& old_calo_shower_step) {

  // Need this if-statement because sometimes the SimParticle that is being Ptr'd to 
  // is not there... The Ptr itself is valid (i.e. old_step.simParticle().isNonnull() returns true)
  // but there is no object there and so when we try to get the id of the SimParticle
  // there is a segfault
  if (old_calo_shower_step.simParticle().get()) {
    CaloShowerStep new_calo_shower_step = old_calo_shower_step;

    _newCaloShowerSteps->push_back(new_calo_shower_step);

    return art::Ptr<mu2e::CaloShowerStep>(_newCaloShowerStepsPID, _newCaloShowerSteps->size()-1, _newCaloShowerStepGetter);
  }
  else {
    return art::Ptr<CaloShowerStep>();
  }
}

void mu2e::FixCaloShowerStepPtrs::copyCaloShowerSim(const mu2e::CaloShowerSim& old_calo_shower_sim, const CaloShowerStepRemap& remap) {

  const auto& caloShowerStepPtrs = old_calo_shower_sim.caloShowerSteps();
  std::vector<art::Ptr<CaloShowerStep> > newCaloShowerStepPtrs;
  for (const auto& i_caloShowerStepPtr : caloShowerStepPtrs) {
    //    std::cout << "AE: " << i_caloShowerStepPtr << std::endl;
    newCaloShowerStepPtrs.push_back(remap.at(i_caloShowerStepPtr));
  }

  CaloShowerSim new_calo_shower_sim = old_calo_shower_sim;
  new_calo_shower_sim.setCaloShowerSteps(newCaloShowerStepPtrs);

  _newCaloShowerSims->push_back(new_calo_shower_sim);
}

void mu2e::FixCaloShowerStepPtrs::copyCaloShowerStepRO(const mu2e::CaloShowerStepRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap) {

  const auto& caloShowerStepPtr = old_calo_shower_step_ro.caloShowerStep();
  CaloShowerStepRO new_calo_shower_step_ro = old_calo_shower_step_ro;

  new_calo_shower_step_ro.setCaloShowerStep(remap.at(caloShowerStepPtr));

  _newCaloShowerStepROs->push_back(new_calo_shower_step_ro);
}

DEFINE_ART_MODULE(mu2e::FixCaloShowerStepPtrs)
