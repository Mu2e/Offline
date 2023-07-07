// Andy Edmonds, 2020
// based on StoppedParticleReactionGun by Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Offline/GeneralUtilities/inc/RSNTIO.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "TH1.h"

namespace mu2e {

  //================================================================
  class MuStopProductsGun : public art::EDProducer {

    typedef RootTreeSampler<IO::StoppedParticleF> RTS;

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter captureProducts{Name("captureProducts"), Comment("Products coincident with captured muon)")};
      fhicl::DelegatedParameter decayProducts{Name("decayProducts"), Comment("Products coincident with decayed muon)")};
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity Level (default = 0)"), 0};
      fhicl::Atom<std::string> material{Name("material"), Comment("Material in which muon is stopped"), "Al"};
      fhicl::Table<RTS::Config> stops{Name("stops"), Comment("Stops ntuple config")};
    };
    typedef art::EDProducer::Table<Config> Parameters;

  private:
    Config conf_;

    int               verbosityLevel_;
    std::string material_;
    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randomFlat_;

    RTS stops_;

    double _decayFraction;
    double _captureFraction;


    std::vector<std::unique_ptr<ParticleGeneratorTool>> _muonDecayGenerators;
    std::vector<std::unique_ptr<ParticleGeneratorTool>> _muonCaptureGenerators;

  public:
    explicit MuStopProductsGun(const Parameters& conf);

    virtual void produce(art::Event& event);
  };

  //================================================================
  MuStopProductsGun::MuStopProductsGun(const Parameters& conf)
    : EDProducer(conf)
    , conf_(conf())
    , verbosityLevel_(conf_.verbosityLevel())
    , material_ (conf_.material())
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randomFlat_(eng_)
    , stops_(eng_, conf_.stops())
    , _decayFraction(GlobalConstantsHandle<PhysicsParams>()->getDecayFraction())
    , _captureFraction(1 - _decayFraction)
  {
    produces<mu2e::GenParticleCollection>();

    if (verbosityLevel_ > 0) {
      std::cout<<"MuStopProductsGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;
      std::cout << "MuStopProductsGun: decayFraction = " << _decayFraction << std::endl;
      std::cout << "MuStopProductsGun: captureFraction = " << _captureFraction << std::endl;
    }

    const auto cap_psets = conf_.captureProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& i_cap_pset : cap_psets) {
      _muonCaptureGenerators.push_back(art::make_tool<ParticleGeneratorTool>(i_cap_pset));
      _muonCaptureGenerators.back()->finishInitialization(eng_, material_);
    }
    const auto decay_psets = conf_.decayProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& i_decay_pset : decay_psets) {
      _muonDecayGenerators.push_back(art::make_tool<ParticleGeneratorTool>(i_decay_pset));
      _muonDecayGenerators.back()->finishInitialization(eng_, material_);
    }
  }

  //================================================================
  void MuStopProductsGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    double rand = randomFlat_.fire();
    if (rand < _decayFraction) {
      for (const auto& i_muonDecayGenerator : _muonDecayGenerators) {
        i_muonDecayGenerator->generate(output, stop);
      }
    }
    else {
      for (const auto& i_muonCaptureGenerator : _muonCaptureGenerators) {
        i_muonCaptureGenerator->generate(output, stop);
      }
    }

    event.put(std::move(output));
  }
  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuStopProductsGun)
