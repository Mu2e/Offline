////////////////////////////////////////////////////////////////////////
// Class:       CondenseRecoTrkCollections
// Plugin Type: prodicer (art v2_06_02)
// File:        CondenseRecoTrkCollections_module.cc
//
// Creates a new RecoTrkBag with new Reco collections that have been reduced
// in size based on a given StrawHitFlag bit (by default onkalseed).
//
// If you want to filter on another StrawHitFlag bit (e.g. you want to store 
// hits near to track hits), then create a new module that will create a
// new StrawHitFlagCollection and the new StrawHitFlag bit. Then you can create a new
// PRecoTrkBag with the new StrawHitFlagCollection and StrawHitFlag bit to this module.
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

#include "DataProducts/inc/ProductBag.hh"

#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"

#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"

namespace mu2e {
  class CondenseRecoTrkCollections;
}


class mu2e::CondenseRecoTrkCollections : public art::EDProducer {
public:
  explicit CondenseRecoTrkCollections(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CondenseRecoTrkCollections(CondenseRecoTrkCollections const &) = delete;
  CondenseRecoTrkCollections(CondenseRecoTrkCollections &&) = delete;
  CondenseRecoTrkCollections & operator = (CondenseRecoTrkCollections const &) = delete;
  CondenseRecoTrkCollections & operator = (CondenseRecoTrkCollections &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void addStrawHitRecoProducts(StrawHitIndex index);

private:

  // the straw hit flag that the new straw hit products need to have (see RecoDataProducts/src/StrawHitFlag.cc for the options)
  std::string _wantedHitFlag;

  // art tags for the input collections
  art::InputTag _trkBagTag;

  // handles to the old collections
  art::Handle<StrawHitFlagCollection> _strawHitFlagsHandle;
  art::Handle<StrawHitCollection> _strawHitsHandle;
  art::Handle<StrawDigiCollection> _strawDigisHandle;
  art::Handle<StrawHitPositionCollection> _strawHitPositionsHandle;

  art::Handle<KalSeedCollection> _kalFinalFitsHandle;
  art::Handle<TrkQualCollection> _trkQualsHandle;
  art::Handle<ProductBag> _trkBagHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<StrawHitFlagCollection> _newStrawHitFlags;
  std::unique_ptr<StrawHitCollection> _newStrawHits;
  std::unique_ptr<StrawDigiCollection> _newStrawDigis;
  std::unique_ptr<StrawHitPositionCollection> _newStrawHitPositions;

  std::unique_ptr<KalSeedCollection> _newKalFinalFits;

  std::map<StrawHitIndex, StrawHitIndex> _oldToNewStrawHitIndexMap;
  int _hitCounter; // keep track of how many straw hit objects we've added to each output collections
};


mu2e::CondenseRecoTrkCollections::CondenseRecoTrkCollections(fhicl::ParameterSet const & pset)
  : _wantedHitFlag(pset.get<std::string>("wantedHitFlag")),
    _trkBagTag(pset.get<art::InputTag>("trkBagTag"))
{
  // Call appropriate produces<>() functions here.
  produces<StrawHitFlagCollection>();
  produces<StrawHitCollection>();
  produces<StrawDigiCollection>();
  produces<StrawHitPositionCollection>();
  
  produces<KalSeedCollection>();
}

void mu2e::CondenseRecoTrkCollections::produce(art::Event & event)
{
  // Implementation of required member function here.

  _newStrawHitFlags = std::unique_ptr<StrawHitFlagCollection>(new StrawHitFlagCollection);
  _newStrawHits = std::unique_ptr<StrawHitCollection>(new StrawHitCollection);
  _newStrawDigis = std::unique_ptr<StrawDigiCollection>(new StrawDigiCollection);
  _newStrawHitPositions = std::unique_ptr<StrawHitPositionCollection>(new StrawHitPositionCollection);
  _newKalFinalFits = std::unique_ptr<KalSeedCollection>(new KalSeedCollection);

  _hitCounter = 0;
  _oldToNewStrawHitIndexMap.clear();

  event.getByLabel(_trkBagTag, _trkBagHandle);
  const auto& trkBag = *_trkBagHandle;
  
  trkBag.getHandle(event, _strawHitFlagsHandle);
  if (!_strawHitFlagsHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find StrawHitFlagCollection in ProductBag\n";
  }

  trkBag.getHandle(event, _strawHitsHandle);
  if (!_strawHitsHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find StrawHitCollection in ProductBag\n";
  }
    
  trkBag.getHandle(event, _strawDigisHandle);
  if (!_strawDigisHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find StrawDigiCollection in ProductBag\n";
  }
  
  trkBag.getHandle(event, _strawHitPositionsHandle);
  if (!_strawHitPositionsHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find StrawHitPositionCollection in ProductBag\n";
  }
  
  trkBag.getHandle(event, _kalFinalFitsHandle);
  if (!_kalFinalFitsHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find KalSeedCollection in ProductBag\n";
  }
  
  trkBag.getHandle(event, _trkQualsHandle);
  if (!_trkQualsHandle.isValid()) {
    throw cet::exception("CondenseRecoTrkCollections") << "Couldn't find TrkQualCollection in ProductBag\n";
  }


  // Loop through the straw hit flag collection
  const auto& strawHitFlags(*_strawHitFlagsHandle);
  for (unsigned int i_straw_hit = 0; i_straw_hit < strawHitFlags.size(); ++i_straw_hit) {
    const mu2e::StrawHitFlag& strawHitFlag = strawHitFlags.at(i_straw_hit);
    //      std::cout << i_straw_hit << ": " << strawHitFlag.stringRep() << std::endl;
    
    // write out the StrawHits that have the StrawHitFlags we want
    mu2e::StrawHitFlag wanted(_wantedHitFlag);
    wanted.merge(StrawHitFlag::onkalseed); // need to have all the hits that are on a KalSeed track (but not necessarily active)
    
    if (strawHitFlag.hasAllProperties(wanted)) {
      StrawHitIndex hit_index = i_straw_hit;
      //      std::cout << i_straw_hit << ": " << strawHitFlag.stringRep() << std::endl;
      addStrawHitRecoProducts(hit_index);
    }
  }
  
  const auto& kalFinalFits(*_kalFinalFitsHandle);
  for (const auto& kalFinalFit : kalFinalFits) {
    
    // Want to check that we have all the hits associated with the KalSeed
    for (const auto& old_hit : kalFinalFit.hits()) {
      try {
	_oldToNewStrawHitIndexMap.at(old_hit.index());
      }
      catch (const std::out_of_range& oor) {
	throw cet::exception("CondenseRecoTrkCollections") << old_hit.index() << " does not exist in _oldToNewStrawHitIndexMap.\n";
      }
    }
    
    KalSeed kal_final_fit(kalFinalFit, _oldToNewStrawHitIndexMap); // copy over the old KalSeed
    _newKalFinalFits->push_back(kal_final_fit);
  }

  // Put everything into the event
  event.put(std::move(_newStrawHitFlags));
  event.put(std::move(_newStrawHits));
  event.put(std::move(_newStrawDigis));
  event.put(std::move(_newStrawHitPositions));

  event.put(std::move(_newKalFinalFits));
}

void mu2e::CondenseRecoTrkCollections::addStrawHitRecoProducts(StrawHitIndex hit_index) {

  //  std::cout << "Will copy out StrawHit with index = " << hit_index << std::endl;
  if (_oldToNewStrawHitIndexMap.find(hit_index) == _oldToNewStrawHitIndexMap.end()) { // only add the straw hit if it's not been seen yet

    mu2e::StrawHitFlag straw_hit_flag = _strawHitFlagsHandle->at(hit_index);
    _newStrawHitFlags->push_back(straw_hit_flag);
  
    mu2e::StrawHit straw_hit = _strawHitsHandle->at(hit_index);
    _newStrawHits->push_back(straw_hit);
  
    mu2e::StrawDigi straw_digi = _strawDigisHandle->at(hit_index);
    _newStrawDigis->push_back(straw_digi);

    mu2e::StrawHitPosition straw_hit_position = _strawHitPositionsHandle->at(hit_index);
    _newStrawHitPositions->push_back(straw_hit_position);

    // Update the map
    _oldToNewStrawHitIndexMap[hit_index] = _hitCounter;
    //    std::cout << hit_index << " --> " << _hitCounter << std::endl;
    ++_hitCounter;
  }
  else {
    // TODO: merge in information from the second StrawHitFlag
  }
}

DEFINE_ART_MODULE(mu2e::CondenseRecoTrkCollections)
