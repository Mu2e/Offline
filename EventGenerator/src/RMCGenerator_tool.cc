// Generator tool to produce radiative muon capture (RMC) events from muon stops
//
// Michael MacKenzie, 2024 (based on DIOGenerator_tool.cc and StoppedMuonRMCGun_module.cc)

#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

namespace mu2e {
  class RMCGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
      fhicl::Atom<std::string> mode{Name("mode"), Comment("Generation mode: \"external\" (real photons), \"internal\" (virtual conversions), or \"physical\" for random sampling given measured rates per muon capture")};
      fhicl::Atom<double> czmin{Name("czmin"), Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double> czmax{Name("czmax"), Comment("Restrict cos(theta_z) maximum"),  1.};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of event kinematics"),  false};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit RMCGenerator(Parameters const& conf) :
      _emass(GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::e_minus).mass()),
      _mumass(GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::mu_minus).mass()),
      _mode(conf().mode()),
      _useRate(_mode == "physical"),
      _external(_useRate || _mode == "external"),
      _czmin(conf().czmin()),
      _czmax(conf().czmax()),
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _internalRate((_external) ? 0. : 1.),
      _makeHistograms(conf().makeHistograms())
    {
      if(!(_mode == "external" || _mode == "internal" || _mode == "physical")) {
        throw cet::exception("BADCONFIG") << "RMCGenerator mode " << _mode << " not defined, only \"external\" \"internal\" or \"physical\" are allowed\n";
      }
      if(_makeHistograms) {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = tfs->mkdir( "RMCGenerator" );
        _hmomentum = tfdir.make<TH1F>("hmomentum", "Produced photon momentum", 60,  0.,  120.  );
        _hCosz     = tfdir.make<TH1F>("hCosz", "Produced photon cos(#theta_{z})", 200,  -1.,  1.  );
        if(!_external || _useRate) { //virtual conversion distributions
          _hEnergyElectron = tfdir.make<TH1F>("hEnergyElectron", "Produced internal conversion electron energy", 60,  0.,  120.  );
          _hEnergyPositron = tfdir.make<TH1F>("hEnergyPositron", "Produced internal conversion positron energy", 60,  0.,  120.  );
          _hMass           = tfdir.make<TH1F>("hMass"          , "M(e+e-) "           , 120,0.,120.);
          _hMassVsE        = tfdir.make<TH2F>("hMassVsE"       , "M(e+e-) vs. E;E (MeV); M(e+e-) (MeV/c^{2})"  , 120,0.,120.,120,0,120);
          _hMassOverE      = tfdir.make<TH1F>("hMassOverE"     , "M(e+e-)/E"          , 100, 0.,1);
          _hy              = tfdir.make<TH1F>("hy"             , "y = (ee-ep)/|pe+pp|", 100,-1.,1.);
        }
      }

      // Validate the input
      const auto s_conf = conf().spectrum.get<fhicl::ParameterSet>();
      const double elow  = s_conf.get<double>("elow");
      const double ehigh = s_conf.get<double>("ehi");
      if(elow < 0.) throw cet::exception("BADCONFIG") << "RMCGenerator minimum energy " << elow
                                                      << " below 0\n";
      if(ehigh > _mumass) throw cet::exception("BADCONFIG") << "RMCGenerator maximum energy " << ehigh
                                                            << " above the muon mass " << _mumass
                                                            << "\n";
      if(ehigh < elow) throw cet::exception("BADCONFIG") << "RMCGenerator minimum energy " << elow
                                                         << " above maximum energy " << ehigh
                                                         << "\n";
      if(_czmin > _czmax || _czmin < -1. || _czmax > 1.) throw cet::exception("BADCONFIG") << "RMCGenerator cos(theta_z) range is not defined\n";

      if(!_external) {
        if(elow < 2.*_emass) throw cet::exception("BADCONFIG") << "RMCGenerator minimum energy " << elow
                                                               << " below the necessary threshold of "
                                                               << 2.*_emass
                                                               << "\n";
      }
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material) override {
      _randomUnitSphereExternal = std::make_unique<RandomUnitSphere>(eng, _czmin, _czmax);
      _randomUnitSphereInternal = std::make_unique<RandomUnitSphere>(eng);
      _randFlat = std::make_unique<CLHEP::RandFlat>(eng);
      _randSpectrum = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
      _muonCaptureSpectrum = std::make_unique<MuonCaptureSpectrum>(_randFlat.get(), _randomUnitSphereInternal.get());
      if(_useRate) {
        _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, GlobalConstantsHandle<PhysicsParams>()->getCaptureRMCRate(material));
        _internalRate = GlobalConstantsHandle<PhysicsParams>()->getCaptureRMCInternalRate(material);
      }
    }

  private:
    const double _emass;
    const double _mumass;

    const std::string _mode; //generation mode: external, internal, or physical
    const bool _useRate; //use the RMC rate to produce N events
    const bool _external; //real or virtual photons if not physical rates
    const double _czmin; //range of cos(theta_z) generated
    const double _czmax;
    BinnedSpectrum _spectrum; //RMC photon spectrum
    std::unique_ptr<MuonCaptureSpectrum> _muonCaptureSpectrum; // internal conversion spectrum
    double _internalRate;

    const bool _makeHistograms;

    std::unique_ptr<RandomUnitSphere>   _randomUnitSphereExternal;
    std::unique_ptr<RandomUnitSphere>   _randomUnitSphereInternal;
    std::unique_ptr<CLHEP::RandFlat>    _randFlat;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
    std::unique_ptr<CLHEP::RandPoissonQ> _randomPoissonQ;

    TH1* _hmomentum;
    TH1* _hCosz;
    TH1* _hEnergyElectron;
    TH1* _hEnergyPositron;
    TH1* _hMass;
    TH2* _hMassVsE;
    TH1* _hMassOverE;
    TH1* _hy; // splitting function
  };


  std::vector<ParticleGeneratorTool::Kinematic> RMCGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    const int n_gen = (_useRate) ? _randomPoissonQ->fire() : 1;
    for(int i_gen = 0; i_gen < n_gen; ++i_gen) {
      // real or virtual photon energy
      const double energy = _spectrum.sample(_randSpectrum->fire());
      if(_makeHistograms) _hmomentum->Fill(energy);
      //determine it's real or virtual conversion
      const bool external = (_useRate) ? _randFlat->fire() > _internalRate : _external;
      if(external) { //real photon spectrum
        const CLHEP::Hep3Vector p3 = _randomUnitSphereExternal->fire(energy);
        const CLHEP::HepLorentzVector fourmom(p3, energy);

        ParticleGeneratorTool::Kinematic k{PDGCode::gamma, ProcessCode::mu2eExternalRMC, fourmom};
        res.emplace_back(k);
        if(_makeHistograms) _hCosz->Fill(fourmom.cosTheta());
      } else { //virtual photon conversion spectrum
        CLHEP::HepLorentzVector p_em, p_ep;
        _muonCaptureSpectrum->getElecPosiVectors(energy, p_em, p_ep);

        ParticleGeneratorTool::Kinematic k_em{PDGCode::e_minus, ProcessCode::mu2eInternalRMC, p_em};
        res.emplace_back(k_em);
        ParticleGeneratorTool::Kinematic k_ep{PDGCode::e_plus , ProcessCode::mu2eInternalRMC, p_ep};
        res.emplace_back(k_ep);
        if(_makeHistograms) { //store virtual kinematics if requested
          _hCosz->Fill((p_em+p_ep).cosTheta());
          _hEnergyElectron->Fill(p_em.e());
          _hEnergyPositron->Fill(p_ep.e());

          const double mass = (p_em+p_ep).m();
          _hMass->Fill(mass);
          _hMassVsE->Fill(energy,mass);
          _hMassOverE->Fill(mass/energy);

          const CLHEP::Hep3Vector p = p_em.vect()+p_ep.vect();
          const double y = (p_em.e()-p_ep.e())/p.mag();
          _hy->Fill(y);
        }
      }
    }

    return res;
  }

  // Legacy implementation of muon stop sampling
  void RMCGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        (d.pdgId == PDGCode::gamma) ? GenId::ExternalRMC : GenId::InternalRMC,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::RMCGenerator)
