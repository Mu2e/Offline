////////////////////////////////////////////////////////////////////////
// Class:       FlagInTimeHits
// Plugin Type: producer (art v2_07_03)
// File:        FlagInTimeHits_module.cc
//
// Flags hits that are on the KalSeed
//
// Generated at Wed Oct 11 15:18:54 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v3_00_01.
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

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

namespace mu2e {
  class FlagInTimeHits;
}


class mu2e::FlagInTimeHits : public art::EDProducer {
public:
  explicit FlagInTimeHits(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlagInTimeHits(FlagInTimeHits const &) = delete;
  FlagInTimeHits(FlagInTimeHits &&) = delete;
  FlagInTimeHits & operator = (FlagInTimeHits const &) = delete;
  FlagInTimeHits & operator = (FlagInTimeHits &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // art tags for the input collections
  art::InputTag _comboHitTag;
  art::InputTag _simParticleTag;

  // Handles to old collections
  art::Handle<ComboHitCollection> _comboHitsHandle;
  art::Handle<SimParticleCollection> _simParticlesHandle;

  // Pointers to new objects that will be put in the event
  std::unique_ptr<ComboHitCollection> _newComboHits;

  double _timeWindow;
  double _mbtime;
  SimParticleTimeOffset _toff;
};


mu2e::FlagInTimeHits::FlagInTimeHits(fhicl::ParameterSet const & pset)
  : _comboHitTag(pset.get<art::InputTag>("comboHitTag")),
    _simParticleTag(pset.get<art::InputTag>("simParticleTag")),
    _timeWindow(pset.get<double>("timeWindow")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  produces<ComboHitCollection>();
}

void mu2e::FlagInTimeHits::produce(art::Event & event)
{
  // Implementation of required member function here.
  _newComboHits = std::unique_ptr<ComboHitCollection>(new ComboHitCollection);
  _toff.updateMap(event);

  // Find the primary particle
  event.getByLabel(_simParticleTag, _simParticlesHandle);
  const auto& simParticles = *_simParticlesHandle;
  double primary_time = 0;
  ConditionsHandle<AcceleratorParams> accPar("ignored");
  _mbtime = accPar->deBuncherPeriod;
  for (const auto& i_sim_particle : simParticles) {
    if (!i_sim_particle.second.isPrimary()) { // don't want to try and dereference an art:Ptr that doesn't exist
      if (i_sim_particle.second.parent()->isPrimary()) { // want an art::Ptr to the primary SimParticle
	primary_time = std::fmod(i_sim_particle.second.startGlobalTime() + _toff.totalTimeOffset(i_sim_particle.second.parent()), _mbtime);
	break;
      }
    }
  }
  if (primary_time != 0) {

    // Now flag the hits that are close in time
    event.getByLabel(_comboHitTag, _comboHitsHandle);
    const auto& comboHits = *_comboHitsHandle;
    for (const auto& i_combo_hit : comboHits) {
      ComboHit new_hit = i_combo_hit;
      
      const double& i_hit_time = new_hit.time();
      if ( std::fabs(i_hit_time - primary_time) <= _timeWindow) {
	new_hit._flag.merge(StrawHitFlag::intime);  // want to keep track of the old straw hit flag
      }
      
      _newComboHits->push_back(new_hit);
    }
  }


  event.put(std::move(_newComboHits));
}

DEFINE_ART_MODULE(mu2e::FlagInTimeHits)
