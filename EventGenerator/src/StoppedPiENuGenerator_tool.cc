///////////////////////////////////////////////////////////////////////////////
// P.Murat: use 'ProcessCode::mu2ePienu' for pi--> e nu decay
//          also use GenId::piEplusNuGun
///////////////////////////////////////////////////////////////////////////////
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
  class StoppedPiEnuGenerator : public ParticleGeneratorTool {

  private:
    PDGCode::type                       _pdgCode;
    double                              _mass;
    BinnedSpectrum                      _spectrum;
    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;

  public:
    struct PhysConfig {
      using Name   = fhicl::Name;
      using Comment= fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("BinnedSpectrum parameters)")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit StoppedPiEnuGenerator(Parameters const& conf) :
      _pdgCode(PDGCode::e_plus),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgCode).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
    {}

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    virtual ProcessCode   processCode() override { return ProcessCode::mu2ePienu; }

    virtual void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      _randSpectrum     = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }
  };

//-----------------------------------------------------------------------------
  std::vector<ParticleGeneratorTool::Kinematic> StoppedPiEnuGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    double e = _spectrum.sample(_randSpectrum->fire());
    double p = sqrt(e*e -_mass*_mass);
    CLHEP::HepLorentzVector fourmom(_randomUnitSphere->fire(p),e);

    ParticleGeneratorTool::Kinematic k{_pdgCode, ProcessCode::mu2ePienu, fourmom};
    res.emplace_back(k);

    return res;
  }

  void StoppedPiEnuGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId, GenId::piEplusNuGun, pos, d.fourmom, stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::StoppedPiEnuGenerator)
