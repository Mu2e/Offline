#include "art/Utilities/ToolMacros.h"

#include "CLHEP/Random/RandFlat.h"

#include "EventGenerator/inc/ParticleGeneratorTool.hh"

#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {
  class MuCap1809keVGammaGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCap1809keVGammaGenerator(Parameters const& conf) :
      _pdgId(PDGCode::gamma),
      _mass(GlobalConstantsHandle<ParticleDataTable>()->particle(_pdgId).ref().mass().value()),
      _genId(GenId::MuCapGammaRayGenTool),
      _energy(GlobalConstantsHandle<PhysicsParams>()->get1809keVGammaEnergy()),
      _intensity(GlobalConstantsHandle<PhysicsParams>()->get1809keVGammaIntensity())
    {

    }
      void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void setEngine(art::RandomNumberGenerator::base_engine_t& eng) {
      _randomUnitSphere = new RandomUnitSphere(eng);
      _randFlat = new CLHEP::RandFlat(eng);
    }

  private:
    PDGCode::type _pdgId;
    double _mass;
    GenId _genId;
    double _energy;
    double _intensity;

    RandomUnitSphere*   _randomUnitSphere;
    CLHEP::RandFlat* _randFlat;
  };

  void MuCap1809keVGammaGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    double rand = _randFlat->fire();
    if (rand < _intensity) {
      double energy = _energy;
      const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      out->emplace_back(_pdgId,
                        _genId,
                        pos,
                        fourmom,
                        stop.t);
    }
  }
}
DEFINE_ART_CLASS_TOOL(mu2e::MuCap1809keVGammaGenerator)
