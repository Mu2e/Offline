#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

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
  class DIOGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<double>       czmin   {Name("czmin")   , Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double>       czmax   {Name("czmax")   , Comment("Restrict cos(theta_z) maximum"),  1.};
      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit DIOGenerator(Parameters const& conf) :
      _pdgId(PDGCode::e_minus),
      _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
      _czmin(conf().czmin()),
      _czmax(conf().czmax()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _flatSpectrum(conf().spectrum.get<fhicl::ParameterSet>().get<std::string>("spectrumShape", "") == "flat")
    {
      if(_czmin > _czmax || _czmin < -1. || _czmax > 1.) throw cet::exception("BADCONFIG") << "DIOGenerator cos(theta_z) range is not defined\n";

      // compute normalization
      double integral(0.0);
      for(size_t ibin=0;ibin < _spectrum.getNbins();++ibin){
        integral += _spectrum.getPDF(ibin);
      }

      auto fullconfig = conf().spectrum.get<fhicl::ParameterSet>();
      _emin = fullconfig.get<double>("elow", _spectrum.getXMin());
      _emax = fullconfig.get<double>("ehi", _spectrum.getXMax());
      _energyFraction = calculateBinnedSpectrumEnergyFraction(fullconfig);

      std::cout << "[" << __func__ << "] Cos(theta_z) min " << _czmin << " max " << _czmax << std::endl;
      std::cout << "[" << __func__ << "] Restricted Spectrum min " << _spectrum.getAbscissa(0) << " max " << _spectrum.getAbscissa(_spectrum.getNbins()-1) << std::endl;
      std::cout << "[" << __func__ << "] Sampled spectrum fraction " << _energyFraction << std::endl;
      std::cout << "[" << __func__ << "] Sampled spectrum fraction (with cos(theta_z)) " << (_energyFraction)*((_czmax - _czmin)/2.) << std::endl;

    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;
    std::unique_ptr<SpectrumConfig> spectrumConfig() override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&, const bool isPrimary) override {
      _isPrimary = isPrimary;
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czmin, _czmax);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
      _randFlat = std::make_unique<CLHEP::RandFlat>(eng);
    }

  private:
    PDGCode::type _pdgId;
    double _mass;

    const double _czmin;
    const double _czmax;
    BinnedSpectrum    _spectrum;
    double _emin;
    double _emax;
    double _energyFraction;
    bool _flatSpectrum;

    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
    std::unique_ptr<CLHEP::RandFlat>    _randFlat;
  };


  std::vector<ParticleGeneratorTool::Kinematic> DIOGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;
    const double r = _energyFraction*(_czmax - _czmin)/2.; // reduce the rate by the spectrum restriction
    if(_isPrimary || _randFlat->fire() <= r) {

      double energy = _spectrum.sample(_randSpectrum->fire());

      const double p = energy * sqrt(1 - std::pow(_mass/energy,2));
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);

      ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2eMuonDecayAtRest, fourmom};
      res.emplace_back(k);
    }

    return res;
  }

  void DIOGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::DIOGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

  std::unique_ptr<SpectrumConfig> DIOGenerator::spectrumConfig() {
    auto config = std::make_unique<SpectrumConfig>();
    config->add_var(SpectrumConfig::RestrictedVar("energy", _energyFraction    , _emin , _emax,
                                                  _flatSpectrum ? SpectrumConfig::Type::kFlat : SpectrumConfig::Type::kPhysical));
    config->add_var(SpectrumConfig::RestrictedVar("cosz"  , (_czmax - _czmin)/2., _czmin, _czmax));
    return config;
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::DIOGenerator)
