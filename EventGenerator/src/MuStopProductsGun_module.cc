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
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "Mu2eUtilities/inc/GenPhysConfig.hh"
#include "EventGenerator/inc/ParticleGeneratorTool.hh"

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
      fhicl::Table<RTS::Config> stops{Name("stops"), Comment("Stops ntuple config")};
    };
    typedef art::EDProducer::Table<Config> Parameters;

  private:
    Config conf_;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };

    struct GenPhysStruct {
      GenPhysStruct(GenPhysConfig conf, art::RandomNumberGenerator::base_engine_t& eng, bool is_poisson_rate = false, double rate = 1, double prob = 1) :
        pdgId(PDGCode::type(conf.pdgId())),
        mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value()),
        spectrumVariable(parseSpectrumVar(conf.spectrumVariable())),
        spectrum(BinnedSpectrum(conf)),
        genId(GenId::findByName(conf.genId())),
        randSpectrum(eng, spectrum.getPDF(), spectrum.getNbins()),
        isPoissonRate(is_poisson_rate),
        emissionRate(rate),
        probability(prob),
        randomFlat(eng),
        randomPoissonQ(eng,rate)
      {    }

      PDGCode::type       pdgId;
      double              mass;
      SpectrumVar       spectrumVariable;
      BinnedSpectrum    spectrum;
      GenId             genId;
      CLHEP::RandGeneral randSpectrum;
      bool isPoissonRate;
      double emissionRate;
      double probability;
      CLHEP::RandFlat randomFlat;
      CLHEP::RandPoissonQ randomPoissonQ;

      int getNGenerated() {
        double rand = randomFlat.fire();
        if (rand < probability) {
          if (isPoissonRate) {
            return randomPoissonQ.fire();
          }
          else {
            return emissionRate;
          }
        }
        else {
          return 0;
        }
      };
    };


    int               verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    RandomUnitSphere   randomUnitSphere_;
    CLHEP::RandFlat randomFlat_;

    RTS stops_;

  private:
    static SpectrumVar    parseSpectrumVar(const std::string& name);
    double                generateEnergy(GenPhysStruct& genPhys);
    void                  generateGenParticles(GenPhysStruct& genPhys, std::unique_ptr<GenParticleCollection>& output);
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
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randomUnitSphere_(eng_)
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
      _muonCaptureGenerators.back()->setEngine(eng_);
    }
    const auto decay_psets = conf_.decayProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& i_decay_pset : decay_psets) {
      _muonDecayGenerators.push_back(art::make_tool<ParticleGeneratorTool>(i_decay_pset));
      _muonDecayGenerators.back()->setEngine(eng_);
    }
  }


  //================================================================
  MuStopProductsGun::SpectrumVar MuStopProductsGun::parseSpectrumVar(const std::string& name) {
    if (name == "totalEnergy"  )  return TOTAL_ENERGY;
    if (name == "kineticEnergy")  return KINETIC_ENERY;
    if (name == "momentum"     )  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"MuStopProductsGun: unknown spectrum variable "<<name<<"\n";
  }


  //================================================================
  void MuStopProductsGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    // // Generate stopped-muon processes here in the future, if needs be
    //    for (auto& i_genPhys : _allStopProducts) {
    //      generateGenParticles(i_genPhys, output);
    //    }
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

//-----------------------------------------------------------------------------
// generate (pseudo-)random particle energy
// the spectrum itself doesn't know whether is stored momentum, kinetic or full
// energy
//-----------------------------------------------------------------------------
  double MuStopProductsGun::generateEnergy(GenPhysStruct& genPhys) {
    double res = genPhys.spectrum.sample(genPhys.randSpectrum.fire());

    if (res < 0.0) {
      throw cet::exception("BADE")<<"MuStopProductsGun: negative energy "<< res <<"\n";
    }

    switch(genPhys.spectrumVariable) {
    case TOTAL_ENERGY  : break;
    case KINETIC_ENERY : res += genPhys.mass; break;
    case MOMENTUM      : res = sqrt(res*res+genPhys.mass*genPhys.mass); break;
    }
    return res;
  }

  // Generate the GenParticles and add them to the output
  void MuStopProductsGun::generateGenParticles(GenPhysStruct& genPhys, std::unique_ptr<GenParticleCollection>& output) {

    const auto& stop = stops_.fire();
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    int n_gen = genPhys.getNGenerated();

    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      const double energy = generateEnergy(genPhys);
      const double p = energy * sqrt(1 - std::pow(genPhys.mass/energy,2));

      CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      output->emplace_back(genPhys.pdgId,
                           genPhys.genId,
                           pos,
                           fourmom,
                           stop.t);

    }
  }
  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuStopProductsGun);
