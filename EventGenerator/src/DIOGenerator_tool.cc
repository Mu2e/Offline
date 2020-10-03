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
  class DIOGenerator : public ParticleGeneratorTool {
  public:
    explicit DIOGenerator(Parameters const& conf) :
      _pdgId(PDGCode::e_minus),
      _mass(GlobalConstantsHandle<ParticleDataTable>()->particle(_pdgId).ref().mass().value()),
      _genId(GenId::DIOGenTool),
      _spectrumVariable(parseSpectrumVar(conf().physics().spectrumVariable())),
      _spectrum(BinnedSpectrum(conf().physics()))
    {

    }
      void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void setEngine(art::RandomNumberGenerator::base_engine_t& eng) {
      _randomUnitSphere = new RandomUnitSphere(eng);
      _randSpectrum = new CLHEP::RandGeneral(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;
    GenId _genId;

    SpectrumVar       _spectrumVariable;
    BinnedSpectrum    _spectrum;

    RandomUnitSphere*   _randomUnitSphere;
    CLHEP::RandGeneral* _randSpectrum;
  };

  void DIOGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    double energy = _spectrum.sample(_randSpectrum->fire());

    switch(_spectrumVariable) {
    case TOTAL_ENERGY  : break;
    case KINETIC_ENERY : energy += _mass; break;
    case MOMENTUM      : energy = sqrt(energy*energy+_mass*_mass); break;
    }

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
DEFINE_ART_CLASS_TOOL(mu2e::DIOGenerator)
