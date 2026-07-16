#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/SpectrumVar.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/DelegatedParameter.h"

namespace mu2e {
  class MuCapNeutronGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
      fhicl::Atom<std::string> spectrumVariable{Name("spectrumVariable"), Comment("Variable that spectrum is defined in (\"totalEnergy\", \"kineticEnergy\", or \"momentum\")")};
      fhicl::Atom<double> czMin{Name("czmin"), Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double> czMax{Name("czmax"), Comment("Restrict cos(theta_z) maximum"),  1.};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCapNeutronGenerator(Parameters const& conf) :
      _pdgId(PDGCode::n0),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _spectrumVariable(parseSpectrumVar(conf().spectrumVariable())),
      _czMin(conf().czMin()),
      _czMax(conf().czMax()),
      _spectrumXMin(_spectrum.getXMin()),
      _spectrumXMax(_spectrum.getXMax()),
      _energyFraction(1.),
      _flatSpectrum(conf().spectrum.get<fhicl::ParameterSet>().get<std::string>("spectrumShape", "") == "flat")
    {
      if(_czMin < -1. || _czMax > 1. || _czMin > _czMax) throw cet::exception("BADCONFIG") << "Cos(theta_z) range unphysical: " << _czMin << " - " << _czMax;
      auto fullconfig = conf().spectrum.get<fhicl::ParameterSet>();
      _energyFraction = calculateBinnedSpectrumEnergyFraction(fullconfig, false);
      std::cout << "[" << __func__ << "] Sampled spectrum fraction " << _energyFraction << std::endl;
      std::cout << "[" << __func__ << "] Sampled spectrum fraction (with cos(theta_z)) " << (_energyFraction)*((_czMax - _czMin)/2.) << std::endl;
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;
    std::unique_ptr<SpectrumConfig> spectrumConfig() override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material, const bool isPrimary) override {
      _isPrimary = isPrimary;
      _rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureNeutronRate(material);
      const double rate = _rate * _energyFraction * (_czMax - _czMin)/2.; // account for potential spectrum restriction in the produced rates
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czMin, _czMax);
      _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, rate);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;
    double _rate = 0.;

    BinnedSpectrum    _spectrum;
    SpectrumVar       _spectrumVariable;
    double _czMin;
    double _czMax;
    double _spectrumXMin;
    double _spectrumXMax;
    double _energyFraction;
    bool _flatSpectrum;

    std::unique_ptr<CLHEP::RandPoissonQ> _randomPoissonQ;
    std::unique_ptr<RandomUnitSphere>    _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral>  _randSpectrum;
  };

  std::vector<ParticleGeneratorTool::Kinematic> MuCapNeutronGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    const int n_gen = (_isPrimary) ? 1 : _randomPoissonQ->fire();
    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      double energy = _spectrum.sample(_randSpectrum->fire());

      switch(_spectrumVariable) {
      case TOTAL_ENERGY  : break;
      case KINETIC_ENERGY : energy += _mass; break;
      case MOMENTUM      : energy = sqrt(energy*energy+_mass*_mass); break;
      }

      const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2eMuonCaptureAtRest, fourmom};
      res.emplace_back(k);
    }

    return res;
  }

  void MuCapNeutronGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapNeutronGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

  std::unique_ptr<SpectrumConfig> MuCapNeutronGenerator::spectrumConfig() {
    auto config = std::make_unique<SpectrumConfig>();

    std::string spectrumVarName("energy");
    switch(_spectrumVariable) {
      case KINETIC_ENERGY: spectrumVarName = "kineticEnergy"; break;
      case MOMENTUM:      spectrumVarName = "momentum";      break;
      default:            spectrumVarName = "energy";        break;
    }

    config->add_var(SpectrumConfig::RestrictedVar(spectrumVarName, _energyFraction, _spectrumXMin, _spectrumXMax,
                                                  _flatSpectrum ? SpectrumConfig::Type::kFlat : SpectrumConfig::Type::kPhysical));
    config->add_var(SpectrumConfig::RestrictedVar("cosz", (_czMax - _czMin)/2., _czMin, _czMax));
    return config;
  }
}
DEFINE_ART_CLASS_TOOL(mu2e::MuCapNeutronGenerator)
