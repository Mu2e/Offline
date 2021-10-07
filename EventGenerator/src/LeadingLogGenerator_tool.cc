#include "art/Utilities/ToolMacros.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "fhiclcpp/types/DelegatedParameter.h"

namespace mu2e {
  class LeadingLogGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit LeadingLogGenerator(Parameters const& conf) :
      _pdgId(conf().pdgId()), 
      _mass(GlobalConstantsHandle<ParticleDataTable>()->particle(_pdgId).ref().mass().value()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
    {
      pid = static_cast<PDGCode::type>(_pdgId);
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material) override {
      _randomUnitSphere = new RandomUnitSphere(eng);
      _randSpectrum = new CLHEP::RandGeneral(eng, _spectrum.getPDF(), _spectrum.getNbins());
      //_endPointEnergy = GlobalConstantsHandle<PhysicsParams>()->GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy(material);
    }

  private:
    int _pdgId;
    ProcessCode _process;
    double _mass;

    BinnedSpectrum    _spectrum;

    RandomUnitSphere*   _randomUnitSphere;
    CLHEP::RandGeneral* _randSpectrum;
    double _endPointEnergy;
    
    PDGCode::type pid;
  };


  std::vector<ParticleGeneratorTool::Kinematic> LeadingLogGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    double energy = _spectrum.sample(_randSpectrum->fire());

    const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
    CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
    CLHEP::HepLorentzVector fourmom(p3, energy);

    if(pid == PDGCode::e_minus){ _process = ProcessCode::mu2eCeMinusLeadingLog; }
    if(pid == PDGCode::e_plus){ _process = ProcessCode::mu2eCePlusLeadingLog; }
    ParticleGeneratorTool::Kinematic k{pid, _process, fourmom};
    res.emplace_back(k);

    return res;
  }

  void LeadingLogGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::CeLeadingLogGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::LeadingLogGenerator)
