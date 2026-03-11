#include "art/Utilities/ToolMacros.h"

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
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCapNeutronGenerator(Parameters const& conf) :
      _pdgId(PDGCode::n0),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _spectrumVariable(parseSpectrumVar(conf().spectrumVariable()))
    {}

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material) override {
      _rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureNeutronRate(material);
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, _rate);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;
    double _rate = 0.;

    BinnedSpectrum    _spectrum;
    SpectrumVar       _spectrumVariable;

    std::unique_ptr<CLHEP::RandPoissonQ> _randomPoissonQ;
    std::unique_ptr<RandomUnitSphere>    _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral>  _randSpectrum;
  };

  std::vector<ParticleGeneratorTool::Kinematic> MuCapNeutronGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    int n_gen = _randomPoissonQ->fire();
    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      double energy = _spectrum.sample(_randSpectrum->fire());

      switch(_spectrumVariable) {
      case TOTAL_ENERGY  : break;
      case KINETIC_ENERY : energy += _mass; break;
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
}
DEFINE_ART_CLASS_TOOL(mu2e::MuCapNeutronGenerator)
