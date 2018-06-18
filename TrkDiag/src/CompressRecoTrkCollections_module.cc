////////////////////////////////////////////////////////////////////////
// Class:       CompressRecoTrkCollections
// Plugin Type: prodicer (art v2_06_02)
// File:        CompressRecoTrkCollections_module.cc
//
// Creates new Reco collections that have been reduced
// in size based on a given StrawHitFlag bit.
//
// If you want to filter on another StrawHitFlag bit (e.g. you want to store 
// hits near to track hits), then create a new module that will create a
// new StrawHitFlagCollection and the new StrawHitFlag bit.
//
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

#include <memory>

#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"

//#include "RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {
  class CompressRecoTrkCollections;
}


class mu2e::CompressRecoTrkCollections : public art::EDProducer {
public:
  explicit CompressRecoTrkCollections(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CompressRecoTrkCollections(CompressRecoTrkCollections const &) = delete;
  CompressRecoTrkCollections(CompressRecoTrkCollections &&) = delete;
  CompressRecoTrkCollections & operator = (CompressRecoTrkCollections const &) = delete;
  CompressRecoTrkCollections & operator = (CompressRecoTrkCollections &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void addStrawHitRecoProducts(StrawHitIndex index);

private:

  // the straw hit flag that the new straw hit products need to have (see RecoDataProducts/src/StrawHitFlag.cc for the options)
  std::string _wantedHitFlag;

  // art tags for the input collections
  art::InputTag _comboHitTag;
  art::InputTag _strawDigiTag;

  // handles to the old collections
  art::Handle<ComboHitCollection> _comboHitsHandle;
  art::Handle<StrawDigiCollection> _strawDigisHandle;

  //  art::Handle<KalSeedCollection> _kalFinalFitsHandle;
  //  art::Handle<TrkQualCollection> _trkQualsHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<ComboHitCollection> _newComboHits;
  std::unique_ptr<StrawDigiCollection> _newStrawDigis;

  //  std::unique_ptr<KalSeedCollection> _newKalFinalFits;

  std::map<StrawHitIndex, StrawHitIndex> _oldToNewStrawHitIndexMap;
  int _hitCounter; // keep track of how many straw hit objects we've added to each output collections
};


mu2e::CompressRecoTrkCollections::CompressRecoTrkCollections(fhicl::ParameterSet const & pset)
  : _wantedHitFlag(pset.get<std::string>("wantedHitFlag")),
    _comboHitTag(pset.get<art::InputTag>("comboHitTag")),
    _strawDigiTag(pset.get<art::InputTag>("strawDigiTag"))
{
  // Call appropriate produces<>() functions here.
  produces<ComboHitCollection>();
  produces<StrawDigiCollection>();
  
  //  produces<KalSeedCollection>();
}

void mu2e::CompressRecoTrkCollections::produce(art::Event & event)
{
  // Implementation of required member function here.

  _newComboHits = std::unique_ptr<ComboHitCollection>(new ComboHitCollection);
  _newStrawDigis = std::unique_ptr<StrawDigiCollection>(new StrawDigiCollection);
  //  _newKalFinalFits = std::unique_ptr<KalSeedCollection>(new KalSeedCollection);

  _hitCounter = 0;
  _oldToNewStrawHitIndexMap.clear();

  event.getByLabel(_strawDigiTag, _strawDigisHandle);

  // Loop through the straw hit flag collection
  event.getByLabel(_comboHitTag, _comboHitsHandle);
  const auto& comboHits(*_comboHitsHandle);
  for (unsigned int i_straw_hit = 0; i_straw_hit < comboHits.size(); ++i_straw_hit) {
    const mu2e::StrawHitFlag& strawHitFlag = comboHits.at(i_straw_hit).flag();
    StrawHitIndex hit_index = i_straw_hit;
    
    // write out the StrawHits that have the StrawHitFlags we want
    if (_wantedHitFlag == "") {
      addStrawHitRecoProducts(hit_index);
    }
    else {
      mu2e::StrawHitFlag wanted(_wantedHitFlag);
      if (strawHitFlag.hasAllProperties(wanted)) {
	addStrawHitRecoProducts(hit_index);
      }
    }
  }
  
  /*  const auto& kalFinalFits(*_kalFinalFitsHandle);
  for (const auto& kalFinalFit : kalFinalFits) {
    
    // Want to check that we have all the hits associated with the KalSeed
    for (const auto& old_hit : kalFinalFit.hits()) {
      try {
	_oldToNewStrawHitIndexMap.at(old_hit.index());
      }
      catch (const std::out_of_range& oor) {
	throw cet::exception("CompressRecoTrkCollections") << old_hit.index() << " does not exist in _oldToNewStrawHitIndexMap.\n";
      }
    }
    
    KalSeed kal_final_fit(kalFinalFit, _oldToNewStrawHitIndexMap); // copy over the old KalSeed
    _newKalFinalFits->push_back(kal_final_fit);
  }
  */

  // Put everything into the event
  event.put(std::move(_newComboHits));
  event.put(std::move(_newStrawDigis));

  //  event.put(std::move(_newKalFinalFits));
}

void mu2e::CompressRecoTrkCollections::addStrawHitRecoProducts(StrawHitIndex hit_index) {

  //  std::cout << "Will copy out StrawHit with index = " << hit_index << std::endl;
  if (_oldToNewStrawHitIndexMap.find(hit_index) == _oldToNewStrawHitIndexMap.end()) { // only add the straw hit if it's not been seen yet

    //    mu2e::StrawHitFlag straw_hit_flag = _strawHitFlagsHandle->at(hit_index);
    //    _newStrawHitFlags->push_back(straw_hit_flag);
  
    mu2e::ComboHit combo_hit = _comboHitsHandle->at(hit_index);
    _newComboHits->push_back(combo_hit);
  
    mu2e::StrawDigi straw_digi = _strawDigisHandle->at(hit_index);
    _newStrawDigis->push_back(straw_digi);

    // Update the map
    _oldToNewStrawHitIndexMap[hit_index] = _hitCounter;
    //    std::cout << hit_index << " --> " << _hitCounter << std::endl;
    ++_hitCounter;
  }
  else {
    // TODO: merge in information from the second StrawHitFlag
  }
}

DEFINE_ART_MODULE(mu2e::CompressRecoTrkCollections)
