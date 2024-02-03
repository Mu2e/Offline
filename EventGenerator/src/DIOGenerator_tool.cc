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
#include "fhiclcpp/types/Atom.h"

namespace mu2e {
  class DIOGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
      fhicl::Atom<double> rad_min{Name("rad_min"),0};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit DIOGenerator(Parameters const& conf) :
      _pdgId(PDGCode::e_minus),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _rad_min(conf().rad_min())
    {}

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randomUnitSphere = new RandomUnitSphere(eng);
      _randSpectrum = new CLHEP::RandGeneral(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;
    BinnedSpectrum    _spectrum;
    double _rad_min;

    RandomUnitSphere*   _randomUnitSphere;
    CLHEP::RandGeneral* _randSpectrum;
  };


  std::vector<ParticleGeneratorTool::Kinematic> DIOGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;
    
    double energy = _spectrum.sample(_randSpectrum->fire());

    const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
    CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
    CLHEP::HepLorentzVector fourmom(p3, energy);

    ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2eMuonDecayAtRest, fourmom};
 
    res.emplace_back(k);
    return res;
  }

  void DIOGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);//position from stopped muon
    const auto daughters = generate();

    for(const auto& d: daughters) {

      
      out->emplace_back(d.pdgId,
                        GenId::DIOGenTool,
                        pos,
                        d.fourmom, //CLHEP::HepLorentzVector type
                        stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::DIOGenerator)
