// Generator tool for apparent baryon number violating muon capture into a light dark sector,
// following the benchmark models of Fox, Hostert, Menzo, Pospelov, and Zupan, arXiv:2407.03450:
//
//   recoil = "chi0"     (aLFC/aBNV):  mu- p --> chi0 chi2,  mu- + N(A,Z) --> chi0 chi2 + N(A-1,Z-1)
//   recoil = "electron" (aLFV/aBNV):  mu- n --> e-   chi2,  mu- + N(A,Z) --> e-   chi2 + N(A-1,Z)
//
// followed by the dark cascade
//
//   chi2 --> chi1 A',  A' --> e+ e-   if m(chi2) > m(chi1) + m(A') (on-shell A')
//   chi2 --> chi1 e+ e-               otherwise (off-shell A', three-body phase space)
//
// chi0 and chi1 are invisible and not written out; the e+ e- pair (and the recoil electron, if
// requested) are written out. Since chi2 is produced nearly at rest, the e+ e- energy sum is
// close to (m2^2 - m1^2 + mA^2)/(2 m2), which can exceed the conversion electron energy.
//
// chi2 is prompt by default. With chi2Lifetime > 0 (proper lifetime, ns), chi2 flies in a straight line
// (it is neutral) from the muon stop and the e+ e- pair is produced at the displaced decay point with
// the corresponding lab-frame delay (Kinematic::positionOffset and timeDelay). Only decays inside a
// fiducial region around the detector (maxTransverseFlight, min/maxLongitudinalFlight, relative to
// the stop) are kept: when generating primaries the whole event is resampled otherwise, and the
// fraction of chi2 decaying inside the region is reported in the SpectrumConfig ("chi2DecayInRegion")
// so samples can be normalized. The A' decay is always prompt.
//
// The mu- + nucleon system follows the configured nuclear model, as in MuCapBNVGenerator (see
// MuonNucleonKinematics.hh and the BNVNuclearModel tables in EventGenerator/fcl/prolog.fcl).
// The masses, nuclear parameters and (for a long-lived chi2) the decay region have no defaults.
// All decays are isotropic in the parent rest frame (no spin correlations or matrix elements).
//
// Michael MacKenzie, 2026

#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"
#include "Offline/EventGenerator/inc/MuonNucleonKinematics.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/TableFragment.h"

#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {
  class MuCapDarkCascadeGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<std::string> recoil{Name("recoil"),
          Comment("Particle produced with chi2: \"chi0\" (mu- p --> chi0 chi2) or \"electron\" (mu- n --> e- chi2)"), "chi0"};
      fhicl::OptionalAtom<double> m0{Name("m0"), Comment("chi0 mass (MeV/c^2), required for (and only allowed with) recoil \"chi0\"")};
      fhicl::Atom<double> m1{Name("m1"), Comment("chi1 mass (MeV/c^2)")};
      fhicl::Atom<double> m2{Name("m2"), Comment("chi2 mass (MeV/c^2)")};
      fhicl::Atom<double> mA{Name("mA"), Comment("Dark photon A' mass (MeV/c^2)")};
      fhicl::TableFragment<MuonNucleonKinematics::NuclearConfig> nuclear;
      fhicl::Atom<double> chi2Lifetime{Name("chi2Lifetime"), Comment("chi2 proper lifetime (ns); 0 for a prompt decay"), 0.};
      fhicl::OptionalAtom<double> maxTransverseFlight{Name("maxTransverseFlight"),
          Comment("Keep chi2 decays within this transverse distance of the muon stop (mm); required if chi2Lifetime > 0")};
      fhicl::OptionalAtom<double> minLongitudinalFlight{Name("minLongitudinalFlight"),
          Comment("Keep chi2 decays at least this far along z from the muon stop (mm; negative is upstream); required if chi2Lifetime > 0")};
      fhicl::OptionalAtom<double> maxLongitudinalFlight{Name("maxLongitudinalFlight"),
          Comment("Keep chi2 decays at most this far along z from the muon stop (mm); required if chi2Lifetime > 0")};
      fhicl::Atom<double> rate{Name("rate"), Comment("Rate per muon capture, only used when not generating as a primary"), 0.};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the generated kinematics"), false};
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCapDarkCascadeGenerator(Parameters const& conf);

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material, const bool isPrimary) override;
    std::unique_ptr<SpectrumConfig> spectrumConfig() override;

  private:
    // Sample the chi2 decay point and lab-frame delay; returns false if it lies outside the fiducial region
    bool sampleDecay(const CLHEP::HepLorentzVector& chi2, CLHEP::Hep3Vector& offset, double& delay);

    bool _electronRecoil;
    double _m0; // mass of the particle produced with chi2 (chi0 or the electron)
    double _m1;
    double _m2;
    double _mA;
    double _me;
    bool _onShell; // chi2 --> chi1 A' is open
    PDGCode::type _nucleonPdg;
    double _nucleonMass;
    MuonNucleonKinematics::NuclearParameters _nuclear;
    double _chi2Lifetime;
    double _maxTransverseFlight = 0.;
    double _minLongitudinalFlight = 0.;
    double _maxLongitudinalFlight = 0.;
    double _rate;
    bool _makeHistograms;
    int _verbosity;

    // Fiducial decay bookkeeping for normalization
    unsigned long _nDecays = 0;
    unsigned long _nDecaysInRegion = 0;

    std::unique_ptr<MuonNucleonKinematics>  _kinematics;
    std::unique_ptr<CLHEP::RandPoissonQ>    _randomPoissonQ;
    std::unique_ptr<CLHEP::RandExponential> _randExp;

    TH1* _hNucleonP = nullptr;
    TH1* _hEChi2 = nullptr;
    TH1* _hEm = nullptr;
    TH1* _hEp = nullptr;
    TH1* _hESum = nullptr;
    TH1* _hMee = nullptr;
    TH1* _hCosOpen = nullptr;
    TH1* _hERecoil = nullptr;
    TH2* _hEpVsEm = nullptr;
    TH1* _hDecayLength = nullptr;
    TH1* _hDecayDelay = nullptr;
    TH2* _hDecayRhoVsZ = nullptr;
  };

  //================================================================
  MuCapDarkCascadeGenerator::MuCapDarkCascadeGenerator(Parameters const& conf) :
    _electronRecoil(conf().recoil() == "electron"),
    _m1(conf().m1()),
    _m2(conf().m2()),
    _mA(conf().mA()),
    _nuclear(MuonNucleonKinematics::parseConfig(conf().nuclear(), "MuCapDarkCascadeGenerator")),
    _chi2Lifetime(conf().chi2Lifetime()),
    _rate(conf().rate()),
    _makeHistograms(conf().makeHistograms()),
    _verbosity(conf().verbosity())
  {
    GlobalConstantsHandle<ParticleDataList> pdt;
    _me = pdt->particle(PDGCode::e_minus).mass();

    if(!(_electronRecoil || conf().recoil() == "chi0")) {
      throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: recoil \"" << conf().recoil()
                                        << "\" not defined, only \"chi0\" or \"electron\" are allowed\n";
    }
    // Charge conservation: mu- p --> chi0 chi2 (neutral), mu- n --> e- chi2
    _nucleonPdg  = (_electronRecoil) ? PDGCode::neutron : PDGCode::proton;
    _nucleonMass = pdt->particle(_nucleonPdg).mass();
    double m0 = 0.;
    const bool hasM0 = conf().m0(m0);
    if(hasM0 == _electronRecoil) {
      throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: m0 must " << ((_electronRecoil) ? "not " : "")
                                        << "be set for recoil \"" << conf().recoil() << "\"\n";
    }
    _m0 = (_electronRecoil) ? _me : m0;

    if(_m0 < 0. || _m1 < 0. || _mA < 0.) throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: negative mass\n";
    if(_m2 <= _m1 + 2.*_me) {
      throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: chi2 --> chi1 e+ e- is closed for m2 = " << _m2
                                        << " and m1 = " << _m1 << "\n";
    }
    if(_mA <= 2.*_me) {
      throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: A' --> e+ e- is closed for mA = " << _mA << "\n";
    }
    _onShell = _m2 > _m1 + _mA;
    if(_chi2Lifetime < 0.) throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: negative chi2 lifetime\n";
    const bool hasRho  = conf().maxTransverseFlight(_maxTransverseFlight);
    const bool hasZMin = conf().minLongitudinalFlight(_minLongitudinalFlight);
    const bool hasZMax = conf().maxLongitudinalFlight(_maxLongitudinalFlight);
    const bool hasRegion = hasRho && hasZMin && hasZMax;
    if(_chi2Lifetime > 0.) {
      if(!hasRegion) {
        throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: a long-lived chi2 requires maxTransverseFlight,"
                                          << " minLongitudinalFlight, and maxLongitudinalFlight\n";
      }
      if(_maxTransverseFlight <= 0. || _minLongitudinalFlight >= _maxLongitudinalFlight) {
        throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: empty chi2 fiducial decay region\n";
      }
    }

    if(_verbosity > 0) {
      std::cout << "[MuCapDarkCascadeGenerator] mu- " << ((_electronRecoil) ? "n --> e- chi2" : "p --> chi0 chi2")
                << ", chi2 --> chi1 " << ((_onShell) ? "(A' --> e+ e-)" : "e+ e- (off-shell A')")
                << " with m0 = " << _m0 << ", m1 = " << _m1 << ", m2 = " << _m2 << ", mA = " << _mA << " MeV" << std::endl;
      if(_onShell) std::cout << "[MuCapDarkCascadeGenerator] Approximate E(e+) + E(e-) for chi2 at rest = "
                             << (_m2*_m2 - _m1*_m1 + _mA*_mA)/(2.*_m2) << " MeV" << std::endl;
      if(_chi2Lifetime > 0.) std::cout << "[MuCapDarkCascadeGenerator] chi2 proper lifetime = " << _chi2Lifetime
                                       << " ns, fiducial decay region: rho < " << _maxTransverseFlight << " mm, "
                                       << _minLongitudinalFlight << " < dz < " << _maxLongitudinalFlight << " mm" << std::endl;
    }

    if(_makeHistograms) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("MuCapDarkCascadeGenerator");
      _hNucleonP = tfdir.make<TH1F>("hNucleonP", "Struck nucleon (residual nucleus recoil) momentum;p (MeV/c)", 100, 0., 400.);
      _hEChi2    = tfdir.make<TH1F>("hEChi2"   , "#chi_{2} kinetic energy;T (MeV)", 100, 0., 100.);
      _hEm       = tfdir.make<TH1F>("hEm"      , "Electron energy (from A');E (MeV)", 200, 0., 200.);
      _hEp       = tfdir.make<TH1F>("hEp"      , "Positron energy;E (MeV)", 200, 0., 200.);
      _hESum     = tfdir.make<TH1F>("hESum"    , "e^{+}e^{-} energy sum;E_{+} + E_{-} (MeV)", 200, 0., 200.);
      _hMee      = tfdir.make<TH1F>("hMee"     , "e^{+}e^{-} invariant mass;M_{ee} (MeV/c^{2})", 200, 0., 200.);
      _hCosOpen  = tfdir.make<TH1F>("hCosOpen" , "e^{+}e^{-} opening angle;cos(#theta_{ee})", 100, -1., 1.);
      _hERecoil  = tfdir.make<TH1F>("hERecoil" , "Recoil (#chi_{0} or e^{-}) energy;E (MeV)", 200, 0., 200.);
      _hEpVsEm   = tfdir.make<TH2F>("hEpVsEm"  , "Positron vs. electron energy;E_{-} (MeV);E_{+} (MeV)", 100, 0., 200., 100, 0., 200.);
      if(_chi2Lifetime > 0.) {
        _hDecayLength = tfdir.make<TH1F>("hDecayLength", "#chi_{2} flight distance (accepted decays);L (mm)", 200, 0., 8000.);
        _hDecayDelay  = tfdir.make<TH1F>("hDecayDelay" , "#chi_{2} decay delay in the lab (accepted decays);#Deltat (ns)", 200, 0., 2000.);
        _hDecayRhoVsZ = tfdir.make<TH2F>("hDecayRhoVsZ", "#chi_{2} decay point relative to the stop;#Deltaz (mm);#rho (mm)",
                                         150, _minLongitudinalFlight, _maxLongitudinalFlight, 70, 0., _maxTransverseFlight);
      }
    }
  }

  //================================================================
  void MuCapDarkCascadeGenerator::finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material, const bool isPrimary) {
    _isPrimary = isPrimary;
    MuonNucleonKinematics::checkMaterial(_nuclear, material, "MuCapDarkCascadeGenerator");
    GlobalConstantsHandle<PhysicsParams> phy;
    const double muonEnergy = phy->getMuonEnergy(material);
    // Nuclear mass of the initial nucleus, from the atomic mass
    const double nucleusMass = phy->getAtomicMass(material)*CLHEP::amu_c2
      - phy->getAtomicNumber(material)*GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::e_minus).mass();
    _kinematics = std::make_unique<MuonNucleonKinematics>(eng, muonEnergy, _nucleonMass, nucleusMass,
                                                          _nuclear);
    _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, _rate);
    _randExp = std::make_unique<CLHEP::RandExponential>(eng);

    const double wmax = _kinematics->maxInvariantMass();
    if(wmax <= _m0 + _m2) {
      throw cet::exception("BADCONFIG") << "MuCapDarkCascadeGenerator: final state mass sum " << _m0 + _m2
                                        << " MeV exceeds the maximum available invariant mass " << wmax << " MeV\n";
    }
    if(_verbosity > 0) {
      std::cout << "[MuCapDarkCascadeGenerator] Bound muon energy = " << muonEnergy
                << " MeV, maximum W = " << wmax << " MeV" << std::endl;
    }
  }

  //================================================================
  std::vector<ParticleGeneratorTool::Kinematic> MuCapDarkCascadeGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic> res;

    const int n_gen = (_isPrimary) ? 1 : _randomPoissonQ->fire();
    for(int i_gen = 0; i_gen < n_gen; ++i_gen) {
      // mu- N --> (chi0 or e-) chi2, resampled if a primary chi2 decays outside the fiducial region
      double pNucleon = 0.;
      std::vector<CLHEP::HepLorentzVector> primary;
      CLHEP::Hep3Vector offset;
      double delay = 0.;
      constexpr int maxTries = 1000000;
      int tries = 0;
      bool inRegion = false;
      while(!inRegion) {
        if(++tries > maxTries) {
          throw cet::exception("RUNTIME") << "MuCapDarkCascadeGenerator: no chi2 decay inside the fiducial region after "
                                          << maxTries << " tries; the chi2 lifetime is too long for this region\n";
        }
        _kinematics->generate({_m0, _m2}, primary, &pNucleon);
        inRegion = sampleDecay(primary[1], offset, delay);
        if(!_isPrimary) break; // pileup-style use: a chi2 decaying outside the region is simply lost
      }
      const CLHEP::HepLorentzVector& recoil = primary[0];
      const CLHEP::HepLorentzVector& chi2   = primary[1];
      if(!inRegion) {
        if(_electronRecoil) res.emplace_back(ParticleGeneratorTool::Kinematic{PDGCode::e_minus, ProcessCode::mu2eMuonBNV, recoil});
        continue;
      }

      // chi2 --> chi1 e+ e-
      CLHEP::HepLorentzVector chi1, em, ep;
      if(_onShell) {
        CLHEP::HepLorentzVector aprime;
        _kinematics->twoBodyDecay(chi2, _m1, _mA, chi1, aprime);
        _kinematics->twoBodyDecay(aprime, _me, _me, em, ep);
      } else {
        std::vector<CLHEP::HepLorentzVector> momenta;
        _kinematics->phaseSpace(chi2, {_m1, _me, _me}, momenta);
        chi1 = momenta[0];
        em   = momenta[1];
        ep   = momenta[2];
      }

      // The recoil is produced at the stop; the e+ e- pair at the chi2 decay point
      if(_electronRecoil) res.emplace_back(ParticleGeneratorTool::Kinematic{PDGCode::e_minus, ProcessCode::mu2eMuonBNV, recoil});
      res.emplace_back(ParticleGeneratorTool::Kinematic{PDGCode::e_minus, ProcessCode::mu2eMuonBNV, em, offset, delay});
      res.emplace_back(ParticleGeneratorTool::Kinematic{PDGCode::e_plus , ProcessCode::mu2eMuonBNV, ep, offset, delay});

      if(_makeHistograms) {
        _hNucleonP->Fill(pNucleon);
        _hEChi2->Fill(chi2.e() - _m2);
        _hEm->Fill(em.e());
        _hEp->Fill(ep.e());
        _hESum->Fill(em.e() + ep.e());
        _hMee->Fill((em + ep).m());
        _hCosOpen->Fill(em.vect().unit().dot(ep.vect().unit()));
        _hERecoil->Fill(recoil.e());
        _hEpVsEm->Fill(em.e(), ep.e());
        if(_chi2Lifetime > 0.) {
          _hDecayLength->Fill(offset.mag());
          _hDecayDelay->Fill(delay);
          _hDecayRhoVsZ->Fill(offset.z(), offset.perp());
        }
      }
    }
    return res;
  }

  //================================================================
  void MuCapDarkCascadeGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapBNVGenTool,
                        pos + d.positionOffset,
                        d.fourmom,
                        stop.t + d.timeDelay);
    }
  }

  //================================================================
  bool MuCapDarkCascadeGenerator::sampleDecay(const CLHEP::HepLorentzVector& chi2, CLHEP::Hep3Vector& offset, double& delay) {
    offset = CLHEP::Hep3Vector(0., 0., 0.);
    delay = 0.;
    if(_chi2Lifetime <= 0.) return true; // prompt

    // Straight-line flight of the neutral chi2: proper time t, displacement (p/m) c t, lab delay (E/m) t
    const double t = _randExp->fire(_chi2Lifetime);
    offset = chi2.vect() * (CLHEP::c_light * t / _m2);
    delay  = chi2.e() * t / _m2;
    const bool inRegion = offset.perp() < _maxTransverseFlight
      && offset.z() > _minLongitudinalFlight && offset.z() < _maxLongitudinalFlight;
    ++_nDecays;
    if(inRegion) ++_nDecaysInRegion;
    return inRegion;
  }

  //================================================================
  std::unique_ptr<SpectrumConfig> MuCapDarkCascadeGenerator::spectrumConfig() {
    auto config = std::make_unique<SpectrumConfig>();
    if(_chi2Lifetime > 0.) {
      // Fraction of chi2 decaying inside the fiducial region: the normalization correction for primaries
      const double fraction = (_nDecays > 0) ? double(_nDecaysInRegion)/_nDecays : 1.;
      config->add_var(SpectrumConfig::RestrictedVar("chi2DecayInRegion", fraction, 0., 1.));
      if(_verbosity > 0) {
        std::cout << "[MuCapDarkCascadeGenerator] chi2 decays inside the fiducial region: " << _nDecaysInRegion
                  << " / " << _nDecays << " = " << fraction << std::endl;
      }
    }
    return config;
  }
}
DEFINE_ART_CLASS_TOOL(mu2e::MuCapDarkCascadeGenerator)
