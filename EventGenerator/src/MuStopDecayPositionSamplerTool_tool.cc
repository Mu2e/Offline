// Sample positions / times of muon decays from presampled muon stops
// Ed Callaghan, 2024

// stl
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// art
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// canvas
#include "canvas/Utilities/InputTag.h"

// cetlib_except
#include "cetlib_except/exception.h"

// clhep
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventGenerator/inc/PositionSamplerTool.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e{
  class MuStopDecayPositionSamplerTool:public PositionSamplerTool{
    public:
      struct Config{
        fhicl::Atom<std::string> atom{
          fhicl::Name("atom"),
          fhicl::Comment("atomic symbol")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      MuStopDecayPositionSamplerTool(const Parameters& config):
          _atom(config().atom()){
        auto handle = GlobalConstantsHandle<PhysicsParams>();
        _lifetime = handle->getDecayTime(_atom);
      };

      void UseRandomEngine(art::RandomNumberGenerator::base_engine_t& engine){
        _uniform = std::make_unique<CLHEP::RandFlat>(engine);
        _exponential = std::make_unique<CLHEP::RandExponential>(engine);
      }

      ParticlePositionPair Sample(const SimParticlePtrVector& particles){
        // ...and then sample one...
        size_t idx = static_cast<size_t>(_uniform->fire() * particles.size());
        // this happens with probability 0, but still
        if (idx == particles.size()){
          idx--;
        }
        auto muon = particles[idx];

        // package position and time together
        auto position = muon->endPosition();
        auto time = muon->endGlobalTime() + this->sample_decay_time();
        auto fourpos = CLHEP::HepLorentzVector(position, time);

        // return as pair
        auto rv = std::make_pair(muon, fourpos);
        return rv;
      };

    protected:
      std::string _atom;
      double _lifetime;
      std::unique_ptr<CLHEP::RandFlat> _uniform;
      std::unique_ptr<CLHEP::RandExponential> _exponential;

      double sample_decay_time(){
        double rv = _exponential->fire(_lifetime);
        return rv;
      };
    private:
      /**/
  };
}; // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::MuStopDecayPositionSamplerTool)
