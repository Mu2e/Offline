// Sophie Middleton, 2021
#include "art/Utilities/ToolMacros.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/DelegatedParameter.h"

#include <iostream>

namespace mu2e {
  class Mu2eXGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit Mu2eXGenerator(Parameters const& conf) :
      _pdgId(PDGCode::e_minus),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _emin(0.),
      _emax(0.),
      _energyFraction(1.),
      _flatSpectrum(conf().spectrum.get<fhicl::ParameterSet>().get<std::string>("spectrumShape", "") == "flat")
    {
      _emin = _spectrum.getXMin();
      _emax = _spectrum.getXMax();

      auto fullconfig = conf().spectrum.get<fhicl::ParameterSet>();
      _emin = fullconfig.get<double>("elow", _spectrum.getXMin());
      _emax = fullconfig.get<double>("ehi", _spectrum.getXMax());
      _energyFraction = calculateBinnedSpectrumEnergyFraction(fullconfig, true); // energy fraction evaluation

      std::cout << "[" << __func__ << "] Restricted Spectrum min " << _spectrum.getAbscissa(0) << " max " << _spectrum.getAbscissa(_spectrum.getNbins()-1) << std::endl;
      std::cout << "[" << __func__ << "] Sampled spectrum fraction " << _energyFraction << std::endl;

    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;
    std::unique_ptr<SpectrumConfig> spectrumConfig() override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&, const bool isPrimary) override {
      _isPrimary = isPrimary;
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;

    BinnedSpectrum    _spectrum;
    double _emin;
    double _emax;
    double _energyFraction;
    bool _flatSpectrum;

    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
  };


  std::vector<ParticleGeneratorTool::Kinematic> Mu2eXGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    double energy = _spectrum.sample(_randSpectrum->fire());

    const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
    CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
    CLHEP::HepLorentzVector fourmom(p3, energy);

    ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::Mu2eX, fourmom};
    res.emplace_back(k);

    return res;
  }

  void Mu2eXGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::Mu2eXGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

  std::unique_ptr<SpectrumConfig> Mu2eXGenerator::spectrumConfig() {
    auto config = std::make_unique<SpectrumConfig>();
    config->add_var(SpectrumConfig::RestrictedVar("energy", _energyFraction, _emin, _emax,
                                                  _flatSpectrum ? SpectrumConfig::Type::kFlat : SpectrumConfig::Type::kPhysical));
    return config;
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::Mu2eXGenerator)
