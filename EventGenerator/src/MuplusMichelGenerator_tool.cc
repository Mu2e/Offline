#include "art/Utilities/ToolMacros.h"
#include <memory>

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

namespace mu2e {
  class MuplusMichelGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuplusMichelGenerator(Parameters const& conf) :
      _pdgId(PDGCode::e_plus),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
    {}

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;

    BinnedSpectrum    _spectrum;

    std::unique_ptr<RandomUnitSphere>  _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
  };


  std::vector<ParticleGeneratorTool::Kinematic> MuplusMichelGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    double energy = _spectrum.sample(_randSpectrum->fire());

    const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
    CLHEP::HepLorentzVector fourmom(_randomUnitSphere->fire(p), energy);

    ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2ePrimary, fourmom};
    res.emplace_back(k);

    return res;
  }

  void MuplusMichelGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuplusMichelGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::MuplusMichelGenerator)
