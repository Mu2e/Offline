// Generates e+/e- from converted photons.
// This generator uses Bethe Heitler model from G4 routines.
// First implementation by M. Mackenzie in SU2020 repository
//
// S. Di Falco, 2023

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Mu2eUtilities/inc/GammaPairConversionSpectrum.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"

// Geant4 includes
#include "Geant4/G4Material.hh"

namespace mu2e {

  //================================================================
  class GammaConvGenerator : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input converted photons.")};
      fhicl::Atom<std::string> stopMaterial{Name("stopMaterial"),Comment("G4 name of the material where the photon has converted (temporary)."),"ST_Wires"};
      fhicl::Atom<unsigned> verbosity{Name("verbosity")};
      fhicl::Atom<bool>   useCorrelatedAngleOverKE{Name("useCorrelatedAngleOverKE"), Comment("Flag to use correlated e+e- cos/KE"), false};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit GammaConvGenerator(const Parameters& conf);

    void produce(art::Event& event) override;
    void beginSubRun(art::SubRun& subrun) override;

    //----------------------------------------------------------------
  private:

    art::ProductToken<SimParticleCollection> const simsToken_;
    std::string stopMaterial_;
    GammaPairConversionSpectrum::materialData materialData_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    bool   useCorrelatedAngleOverKE_; //correlate the cos/ke for the e+ and e-

    GammaPairConversionSpectrum spectrum_; //pair production spectrum
    enum {kMaxConversionMaterialElements = 10};
    ProcessCode process;
  };

  //================================================================
  GammaConvGenerator::GammaConvGenerator(const Parameters& conf)
    : EDProducer{conf}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , stopMaterial_{conf().stopMaterial()}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , useCorrelatedAngleOverKE_(conf().useCorrelatedAngleOverKE())
    , spectrum_                (&randFlat_, useCorrelatedAngleOverKE_)
  {
    produces<mu2e::StageParticleCollection>();
    process = ProcessCode::mu2eGammaConversion;
  }

  //================================================================
  void GammaConvGenerator::beginSubRun( art::SubRun &subrun){

    G4String const& stopmat=stopMaterial_;
    G4Material*stopMat = findMaterialOrThrow(stopmat);
    auto elements= stopMat->GetElementVector();
    auto const& eleFracs= stopMat->GetFractionVector();
    size_t ele_size=elements->size();
    if (ele_size > kMaxConversionMaterialElements) {
      throw cet::exception("TOO MANY MATERIAL ELEMENTS")
       << ele_size << " while maximum allowed is " << kMaxConversionMaterialElements << "\n";
    }
    //create a corresponding material
    for (size_t iele=0;iele<ele_size;++iele){
      materialData_.elementDatas.push_back(spectrum_.GetElementMap()[(int)elements->at(iele)->GetZ()]);
      materialData_.elementFractions.push_back(eleFracs[iele]);
      if (verbosity_>1){
        std::cout << "Z= " << elements->at(iele)->GetZ() << " frac=" << eleFracs[iele] << std::endl;
      }
    }

  }
  //================================================================
  void GammaConvGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto gammas = simParticleList(simh,PDGCode::gamma,ProcessCode::conv);

    if(gammas.empty()) {
      throw   cet::exception("NO GAMMAS")
        <<"GammaConvGenerator::produce(): no suitable converted gamma in the input SimParticleCollection\n";
    }

    // Take the converted gamma as SimParticle
    //
    const auto gammapart = gammas.at(randFlat_.fireInt((long)gammas.size()));
    const auto p4_gamma = gammapart->endMomentum();
    const auto t_gammastop = gammapart->endGlobalTime();
    if (verbosity_>2){
      std::cout << "GAMMA 4-MOMENTUM: " << p4_gamma << std::endl;
      std::cout << "END VOLUME INDEX: " << gammapart->endVolumeIndex() << std::endl;
    }
    CLHEP::HepLorentzVector p4_eminus, p4_eplus;
    //sample the spectrum
    spectrum_.fire(p4_gamma, materialData_, p4_eminus, p4_eplus);

    output->emplace_back(gammapart,
                         process,
                         PDGCode::e_minus,
                         gammapart->endPosition(),
                         p4_eminus,
                         t_gammastop
                         );
    output->emplace_back(gammapart,
                         process,
                         PDGCode::e_plus,
                         gammapart->endPosition(),
                         p4_eplus,
                         t_gammastop
                         );

    if (verbosity_>1){
      std::cout << "p(gamma)=(" << p4_gamma.vect().x() << "," << p4_gamma.vect().y() << "," << p4_gamma.vect().z() << "," << p4_gamma.e() << ")" << std::endl;
      std::cout << "p(eminus)=(" << p4_eminus.vect().x() << "," << p4_eminus.vect().y() << "," << p4_eminus.vect().z() << "," << p4_eminus.e() << ")" << std::endl;
      std::cout << "p(eplus)=(" << p4_eplus.vect().x() << "," << p4_eplus.vect().y() << "," << p4_eplus.vect().z() << "," << p4_eplus.e() << ")" << std::endl;
    }

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConvGenerator)
