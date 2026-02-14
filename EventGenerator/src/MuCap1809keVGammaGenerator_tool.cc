#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandFlat.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {
  class MuCap1809keVGammaGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> czMin  {Name("czmin")  , Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double> czMax  {Name("czmax")  , Comment("Restrict cos(theta_z) maximum"),  1.};
      fhicl::Atom<bool>   fireAll{Name("fireAll"), Comment("Add a photon to all events, otherwise use the branching fraction"),  false};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCap1809keVGammaGenerator(Parameters const& conf) :
      _czMin(conf().czMin()),
      _czMax(conf().czMax()),
      _fireAll(conf().fireAll()),
      _pdgId(PDGCode::gamma),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _randomUnitSphere(nullptr),
      _randFlat(nullptr)
    {
      if(_czMin < -1. || _czMax > 1. || _czMin > _czMax) throw cet::exception("BADCONFIG") << "Cos(theta_z) range unphysical: " << _czMin << " - " << _czMax;
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material) override {
      _energy = GlobalConstantsHandle<PhysicsParams>()->get1809keVGammaEnergy(material);
      _intensity = GlobalConstantsHandle<PhysicsParams>()->get1809keVGammaIntensity(material);
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czMin, _czMax);
      _randFlat = std::make_unique<CLHEP::RandFlat>(eng);
    }

  private:
    double _czMin;
    double _czMax;
    bool   _fireAll;
    PDGCode::type _pdgId;
    double _mass;
    double _energy = 0.;
    double _intensity = 0.;

    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandFlat> _randFlat;
  };

  std::vector<ParticleGeneratorTool::Kinematic> MuCap1809keVGammaGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    if (_fireAll || _randFlat->fire() < _intensity) {
      const double momentum = _energy * sqrt(1 - std::pow(_mass/_energy,2));
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(momentum);
      CLHEP::HepLorentzVector fourmom(p3, _energy);
      ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2eMuonCaptureAtRest, fourmom};
      res.emplace_back(k);
    }

    return res;
  }

  void MuCap1809keVGammaGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapGammaRayGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::MuCap1809keVGammaGenerator)
