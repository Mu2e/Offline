#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/DelegatedParameter.h"

namespace mu2e {
  class MuCapPhotonGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
      fhicl::Atom<double> czMin{Name("czmin"), Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double> czMax{Name("czmax"), Comment("Restrict cos(theta_z) maximum"),  1.};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCapPhotonGenerator(Parameters const& conf) :
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _czMin(conf().czMin()),
      _czMax(conf().czMax())
    {
      if(_czMin < -1. || _czMax > 1. || _czMin > _czMax) throw cet::exception("BADCONFIG") << "Cos(theta_z) range unphysical: " << _czMin << " - " << _czMax;
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material) override {
      _rate = GlobalConstantsHandle<PhysicsParams>()->getCapturePhotonRate(material);
      const double rate = _rate * (_czMax - _czMin)/2.; // accouunt for potential cz selection in the produced rates
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czMin, _czMax);
      _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, rate);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    double _rate = 0.;

    BinnedSpectrum    _spectrum;
    double _czMin;
    double _czMax;

    std::unique_ptr<CLHEP::RandPoissonQ> _randomPoissonQ;
    std::unique_ptr<RandomUnitSphere>    _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral>  _randSpectrum;
  };


  std::vector<ParticleGeneratorTool::Kinematic> MuCapPhotonGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    int n_gen = _randomPoissonQ->fire();
    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      double energy = _spectrum.sample(_randSpectrum->fire());
      const double p = energy;
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      ParticleGeneratorTool::Kinematic k{PDGCode::gamma, ProcessCode::mu2eMuonCaptureAtRest, fourmom};
      res.emplace_back(k);
    }

    return res;
  }

  void MuCapPhotonGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapPhotonGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }
}
DEFINE_ART_CLASS_TOOL(mu2e::MuCapPhotonGenerator)
