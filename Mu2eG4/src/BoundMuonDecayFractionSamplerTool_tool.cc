// Sample elements from a material based on the relative rates of muons decaying
// after capturing atomically. Each weight is the product of the atomic volume
// and the bound-muon decay fraction.
// Ed Callaghan, 2024

// stl
#include <memory>
#include <string>

// art
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/make_tool.h"

// mu2e
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eG4/inc/AtomicVolumeSamplerTool.hh"
#include "Offline/Mu2eG4/inc/WeightedElementSamplerTool.hh"

namespace mu2e{
  class BoundMuonDecayFractionSamplerTool: public WeightedElementSamplerTool{
    public:
      struct Config{
        fhicl::Atom<std::string> name{
          fhicl::Name("material"),
          fhicl::Comment("Material name")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      BoundMuonDecayFractionSamplerTool(const Parameters&);
     ~BoundMuonDecayFractionSamplerTool(){ /**/ };
      virtual void UseRandomEngine(art::RandomNumberGenerator::base_engine_t&);
    protected:
      std::unique_ptr<AtomicVolumeSamplerTool> _volume_sampler;

      void finish_initialize();
    private:
      /**/
  };

  BoundMuonDecayFractionSamplerTool::BoundMuonDecayFractionSamplerTool(const Parameters& config): WeightedElementSamplerTool(config().name()){
    // initialize atomic volume sampler using same (simple) configuration
    auto pset = config.get_PSet();
    pset.put_or_replace("tool_type", "AtomicVolumeSamplerTool");
    _volume_sampler = art::make_tool<AtomicVolumeSamplerTool>(pset);
  }

  void BoundMuonDecayFractionSamplerTool::finish_initialize(){
    // first, propagate to atomic volume sampler, to calculate capture weights
    // must call a sample here to fake the sampler into post-initializing
    _volume_sampler->Sample();

    // next, adjust the weights by the asymmetric decay fractions
    auto handle = GlobalConstantsHandle<PhysicsParams>();
    auto material = G4Material::GetMaterial(_name, false);
    auto elements = *material->GetElementVector();
    for (size_t i = 0 ; i < elements.size() ; i++){
      auto element = elements[i];
      std::string name = static_cast<std::string>(element->GetName());
      double fraction = handle->getDecayFraction(name);
      _weights[name] = _volume_sampler->GetWeight(name) * fraction;
    }

    // finally, re-normalize the weight vector
    this->normalize_weights();
  }

  void BoundMuonDecayFractionSamplerTool::UseRandomEngine(art::RandomNumberGenerator::base_engine_t& engine){
    // propagate engine to volume sampler
    _volume_sampler->UseRandomEngine(engine);
    // initialize our own distributions
    this->WeightedElementSamplerTool::UseRandomEngine(engine);
  }
}

DEFINE_ART_CLASS_TOOL(mu2e::BoundMuonDecayFractionSamplerTool)
