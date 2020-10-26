#include "art/Utilities/ToolMacros.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "EventGenerator/inc/ParticleGeneratorTool.hh"

#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {
  class MuCapPhotonGenerator : public ParticleGeneratorTool {
  public:
    explicit MuCapPhotonGenerator(Parameters const& conf) :// fhicl::ParameterSet const& pset)
      _pdgId(PDGCode::gamma),
      _genId(GenId::MuCapPhotonGenTool),
      _rate(GlobalConstantsHandle<PhysicsParams>()->getCapturePhotonRate()),
      _spectrum(BinnedSpectrum(conf().physics()))
    {

    }
      void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void setEngine(art::RandomNumberGenerator::base_engine_t& eng) {
      _randomUnitSphere = new RandomUnitSphere(eng);
      _randomPoissonQ = new CLHEP::RandPoissonQ(eng, _rate);
      _randSpectrum = new CLHEP::RandGeneral(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    GenId _genId;
    double _rate;

    BinnedSpectrum    _spectrum;


    CLHEP::RandPoissonQ* _randomPoissonQ;
    RandomUnitSphere*   _randomUnitSphere;
    CLHEP::RandGeneral* _randSpectrum;
  };

  void MuCapPhotonGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    int n_gen = _randomPoissonQ->fire();
    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      double energy = _spectrum.sample(_randSpectrum->fire());

      const double p = energy;
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
DEFINE_ART_CLASS_TOOL(mu2e::MuCapPhotonGenerator)
