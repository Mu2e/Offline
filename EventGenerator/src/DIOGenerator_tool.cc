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
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
    {
      if(_czmin > _czmax || _czmin < -1. || _czmax > 1.) throw cet::exception("BADCONFIG") << "DIOGenerator cos(theta_z) range is not defined\n";

      // compute normalization
      double integral(0.0);
      for(size_t ibin=0;ibin < _spectrum.getNbins();++ibin){
        integral += _spectrum.getPDF(ibin);
      }

      auto fullconfig = conf().spectrum.get<fhicl::ParameterSet>();
      fullconfig.erase(std::string("elow"));
      fullconfig.erase(std::string("ehi"));
      fullconfig.put(std::string("elow"),double(0.0));
      fullconfig.put(std::string("ehi"),double(0.0));
      BinnedSpectrum fullspect(fullconfig);
      double fullintegral(0.0);
      for(size_t ibin=0;ibin < fullspect.getNbins();++ibin){
        fullintegral += fullspect.getPDF(ibin);
      }
      // correct for the missing prediction near threshold, assuming the rate falls to 0 linearly.
      double pmin = _spectrum.getAbscissa(0);
      double pdfmin = _spectrum.getPDF(0);
      double binsize = _spectrum.getBinWidth();
      fullintegral += 0.5*pdfmin*pmin/binsize;
      std::cout << "Cos(theta_z) min " << _czmin << " max " << _czmax << std::endl;
      std::cout << "Restricted Spectrum min " << _spectrum.getAbscissa(0) << " max " << _spectrum.getAbscissa(_spectrum.getNbins()-1) << std::endl;
      std::cout << "Full Spectrum min " << fullspect.getAbscissa(0) << " max " << fullspect.getAbscissa(fullspect.getNbins()-1) << std::endl;
      std::cout << "Restricted Spectrum integral " << integral << std::endl;
      std::cout << "Restricted Spectrum integral*cos(theta_z) restriction " << integral*((_czmax - _czmin)/2.) << std::endl;
      std::cout << "Full Spectrum integral " << fullintegral << std::endl;
      std::cout << "Sampled spectrum fraction " << integral/fullintegral << std::endl;
      std::cout << "Sampled spectrum fraction (with cos(theta_z)) " << (integral/fullintegral)*((_czmax - _czmin)/2.) << std::endl;

    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czmin, _czmax);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }

  private:
    PDGCode::type _pdgId;
    double _mass;

    const double _czmin;
    const double _czmax;
    BinnedSpectrum    _spectrum;

    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
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

}
DEFINE_ART_CLASS_TOOL(mu2e::DIOGenerator)
