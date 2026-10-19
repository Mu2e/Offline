// Generator tool for baryon number violating (BNV) muon-nucleon annihilation in a muonic atom:
//
//   mu- + N(A,Z) --> X_1 + ... + X_n + N(A-1,Z-1)   (bound proton annihilated)
//   mu- + N(A,Z) --> X_1 + ... + X_n + N(A-1,Z)     (bound neutron annihilated)
//
// The final state X_1...X_n is configurable: visible daughters are given by PDG id and are
// written out, invisible daughters (e.g. a dark fermion, a dark scalar, or a neutrino) are given
// by mass only and enter the kinematics but are not written out.
//
// Kinematics: the muon is taken at rest with its bound 1s total energy (PhysicsParams muonEnergy)
// and the struck nucleon follows the configured nuclear model (no defaults, see the
// BNVNuclearModel tables in EventGenerator/fcl/prolog.fcl). The final state is distributed
// according to pure n-body Lorentz-invariant phase space at the invariant mass W of the
// annihilating system (see MuonNucleonKinematics.hh).
// No matrix element and no intra-nuclear final-state interactions of the produced hadrons are
// applied. For the light dark sector final states see Fox et al., arXiv:2407.03450, and for the
// B = L = 1 boson see Heeck, arXiv:2009.01256.
//
// Michael MacKenzie, 2026 (based on MuCapProtonGenerator_tool.cc and RMCGenerator_tool.cc)

#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"

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
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableFragment.h"

#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {
  class MuCapBNVGenerator : public ParticleGeneratorTool {
  public:
    struct InvisibleConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> mass{Name("mass"), Comment("Invisible particle mass (MeV/c^2)")};
      fhicl::Atom<bool> fermion{Name("fermion"), Comment("Is the invisible particle a fermion (used in the spin-statistics check)")};
    };

    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<int> daughters{Name("daughters"), Comment("PDG ids of the visible final state particles")};
      fhicl::OptionalSequence<fhicl::Table<InvisibleConfig>> invisibles{Name("invisibles"),
          Comment("Invisible final state particles, included in the kinematics but not written out")};
      fhicl::Atom<std::string> nucleon{Name("nucleon"), Comment("Annihilated bound nucleon: \"proton\" or \"neutron\""), "proton"};
      fhicl::TableFragment<MuonNucleonKinematics::NuclearConfig> nuclear;
      fhicl::Atom<double> rate{Name("rate"), Comment("BNV rate per muon capture, only used when not generating as a primary"), 0.};
      fhicl::Atom<bool> requireSpinStatistics{Name("requireSpinStatistics"),
          Comment("Throw if the final state has an odd number of fermions (angular momentum cannot be conserved)"), true};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the generated kinematics"), false};
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuCapBNVGenerator(Parameters const& conf);

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material, const bool isPrimary) override;

  private:
    struct Daughter {
      PDGCode::type pdgId;
      double mass;
      bool visible;
    };

    static bool isFermion(int pdgId);
    static int baryonNumber(int pdgId);
    static int leptonNumber(int pdgId);

    std::vector<Daughter> _daughters;
    std::vector<double> _masses;
    PDGCode::type _nucleonPdg;
    double _nucleonMass;
    MuonNucleonKinematics::NuclearParameters _nuclear;
    double _rate;
    bool _makeHistograms;
    int _verbosity;

    double _massSum = 0.;

    std::unique_ptr<MuonNucleonKinematics> _kinematics;
    std::unique_ptr<CLHEP::RandPoissonQ>   _randomPoissonQ;

    std::vector<TH1*> _hDaughterP;
    std::vector<TH1*> _hDaughterCosz;
    TH1* _hNucleonP = nullptr;
    TH1* _hW = nullptr;
    TH1* _hEVis = nullptr;
    TH1* _hPVis = nullptr;
    TH1* _hMVis = nullptr;
    TH2* _hP0VsP1 = nullptr;
  };

  //================================================================
  MuCapBNVGenerator::MuCapBNVGenerator(Parameters const& conf) :
    _nuclear(MuonNucleonKinematics::parseConfig(conf().nuclear(), "MuCapBNVGenerator")),
    _rate(conf().rate()),
    _makeHistograms(conf().makeHistograms()),
    _verbosity(conf().verbosity())
  {
    GlobalConstantsHandle<ParticleDataList> pdt;

    // The annihilated nucleon
    if(conf().nucleon() == "proton")       _nucleonPdg = PDGCode::proton;
    else if(conf().nucleon() == "neutron") _nucleonPdg = PDGCode::neutron;
    else throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: nucleon \"" << conf().nucleon()
                                           << "\" not defined, only \"proton\" or \"neutron\" are allowed\n";
    _nucleonMass = pdt->particle(_nucleonPdg).mass();


    // The final state particles
    int nFermions = 0;
    double charge = 0.;
    int dB = -baryonNumber(_nucleonPdg);
    int dL = -leptonNumber(PDGCode::mu_minus);
    for(const int pdg : conf().daughters()) {
      const auto& pd = pdt->particle(pdg);
      _daughters.push_back(Daughter{static_cast<PDGCode::type>(pdg), pd.mass(), true});
      nFermions += isFermion(pdg);
      charge += pd.charge();
      dB += baryonNumber(pdg);
      dL += leptonNumber(pdg);
    }
    std::vector<InvisibleConfig> invisibles;
    if(conf().invisibles(invisibles)) {
      for(const auto& inv : invisibles) {
        if(inv.mass() < 0.) throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: negative invisible mass\n";
        _daughters.push_back(Daughter{PDGCode::unknown, inv.mass(), false});
        nFermions += inv.fermion();
      }
    }
    if(_daughters.size() < 2) {
      throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: at least two final state particles are required\n";
    }

    // Conservation checks: the invisibles are taken to be neutral
    const double initialCharge = pdt->particle(PDGCode::mu_minus).charge() + pdt->particle(_nucleonPdg).charge();
    if(std::abs(charge - initialCharge) > 1.e-3) {
      throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: final state charge " << charge
                                        << " does not match the initial state charge " << initialCharge << "\n";
    }
    // mu- + nucleon is two fermions, so the final state (the residual nucleus flips statistics along
    // with the removed nucleon, so it cancels) must contain an even number of fermions
    const bool spinStatisticsOK = (nFermions % 2) == 0;
    if(!spinStatisticsOK) {
      if(conf().requireSpinStatistics()) {
        throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: final state has " << nFermions
                                          << " fermion(s); mu- + nucleon annihilation requires an even number."
                                          << " Add an invisible fermion (e.g. a neutrino) or set requireSpinStatistics: false"
                                          << " to generate this final state as a kinematic benchmark only\n";
      }
      std::cout << "[MuCapBNVGenerator] WARNING: final state has an odd number of fermions (" << nFermions
                << ") and violates angular momentum conservation; use only as a kinematic benchmark\n";
    }

    _massSum = 0.;
    for(const auto& d : _daughters) {
      _masses.push_back(d.mass);
      _massSum += d.mass;
    }

    if(_verbosity > 0) {
      std::cout << "[MuCapBNVGenerator] mu- + " << conf().nucleon() << " --> ";
      for(const auto& d : _daughters) {
        if(d.visible) std::cout << pdt->particle(d.pdgId).name() << " ";
        else          std::cout << "invisible(m=" << d.mass << ") ";
      }
      std::cout << "\n[MuCapBNVGenerator] Visible Delta B = " << dB << ", visible Delta L = " << dL
                << " (not counting invisibles), final state fermions = " << nFermions << std::endl;
      std::cout << "[MuCapBNVGenerator] Nuclear model = " << conf().nuclear().nuclearModel() << " for " << _nuclear.material
                << " (k_F = " << _nuclear.fermiMomentum << " MeV/c, removal energy = " << _nuclear.removalEnergy << " MeV)" << std::endl;
    }

    if(_makeHistograms) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("MuCapBNVGenerator");
      for(size_t i = 0; i < _daughters.size(); ++i) {
        if(!_daughters[i].visible) continue;
        const std::string name = pdt->particle(_daughters[i].pdgId).name();
        _hDaughterP   .push_back(tfdir.make<TH1F>(Form("hP_%zu", i), Form("Daughter %zu (%s) momentum;p (MeV/c)", i, name.c_str()), 240, 0., 600.));
        _hDaughterCosz.push_back(tfdir.make<TH1F>(Form("hCosz_%zu", i), Form("Daughter %zu (%s) cos(#theta_{z});cos(#theta_{z})", i, name.c_str()), 100, -1., 1.));
      }
      _hNucleonP = tfdir.make<TH1F>("hNucleonP", "Struck nucleon (residual nucleus recoil) momentum;p (MeV/c)", 100, 0., 400.);
      _hW        = tfdir.make<TH1F>("hW"       , "Invariant mass of all annihilation products;W (MeV/c^{2})", 200, 900., 1100.);
      _hEVis     = tfdir.make<TH1F>("hEVis"    , "Total visible energy;E_{vis} (MeV)", 240, 0., 1200.);
      _hPVis     = tfdir.make<TH1F>("hPVis"    , "Total visible momentum;p_{vis} (MeV/c)", 200, 0., 1000.);
      _hMVis     = tfdir.make<TH1F>("hMVis"    , "Visible invariant mass;M_{vis} (MeV/c^{2})", 240, 0., 1200.);
      _hP0VsP1   = tfdir.make<TH2F>("hP0VsP1"  , "Daughter 0 vs. daughter 1 momentum;p_{1} (MeV/c);p_{0} (MeV/c)", 120, 0., 600., 120, 0., 600.);
    }
  }

  //================================================================
  void MuCapBNVGenerator::finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& material, const bool isPrimary) {
    _isPrimary = isPrimary;
    MuonNucleonKinematics::checkMaterial(_nuclear, material, "MuCapBNVGenerator");
    GlobalConstantsHandle<PhysicsParams> phy;
    const double muonEnergy = phy->getMuonEnergy(material); // bound muon total energy
    // Nuclear mass of the initial nucleus, from the atomic mass
    const double nucleusMass = phy->getAtomicMass(material)*CLHEP::amu_c2
      - phy->getAtomicNumber(material)*GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::e_minus).mass();
    _kinematics = std::make_unique<MuonNucleonKinematics>(eng, muonEnergy, _nucleonMass, nucleusMass,
                                                          _nuclear);
    _randomPoissonQ = std::make_unique<CLHEP::RandPoissonQ>(eng, _rate);

    const double wmax = _kinematics->maxInvariantMass();
    if(wmax <= _massSum) {
      throw cet::exception("BADCONFIG") << "MuCapBNVGenerator: final state mass sum " << _massSum
                                        << " MeV exceeds the maximum available invariant mass " << wmax << " MeV\n";
    }
    if(_verbosity > 0) {
      std::cout << "[MuCapBNVGenerator] Bound muon energy = " << muonEnergy
                << " MeV, maximum W = " << wmax << " MeV" << std::endl;
    }
  }

  //================================================================
  std::vector<ParticleGeneratorTool::Kinematic> MuCapBNVGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic> res;

    const int n_gen = (_isPrimary) ? 1 : _randomPoissonQ->fire();
    for(int i_gen = 0; i_gen < n_gen; ++i_gen) {
      double pNucleon = 0.;
      std::vector<CLHEP::HepLorentzVector> momenta;
      _kinematics->generate(_masses, momenta, &pNucleon);

      CLHEP::HepLorentzVector visible, total;
      for(size_t i = 0; i < _daughters.size(); ++i) {
        total += momenta[i];
        if(!_daughters[i].visible) continue;
        visible += momenta[i];
        res.emplace_back(ParticleGeneratorTool::Kinematic{_daughters[i].pdgId, ProcessCode::mu2eMuonBNV, momenta[i]});
      }

      if(_makeHistograms) {
        size_t ivis = 0;
        for(size_t i = 0; i < _daughters.size(); ++i) {
          if(!_daughters[i].visible) continue;
          _hDaughterP[ivis]->Fill(momenta[i].vect().mag());
          _hDaughterCosz[ivis]->Fill(momenta[i].vect().cosTheta());
          ++ivis;
        }
        _hNucleonP->Fill(pNucleon);
        _hW->Fill(total.m());
        _hEVis->Fill(visible.e());
        _hPVis->Fill(visible.vect().mag());
        _hMVis->Fill(visible.m());
        _hP0VsP1->Fill(momenta[1].vect().mag(), momenta[0].vect().mag());
      }
    }
    return res;
  }

  //================================================================
  void MuCapBNVGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapBNVGenTool,
                        pos,
                        d.fourmom,
                        stop.t);
    }
  }

  //================================================================
  bool MuCapBNVGenerator::isFermion(int pdgId) {
    const int id = std::abs(pdgId);
    if(id >= 11 && id <= 18) return true;           // leptons
    if(id >= 1000000000) return (id/10) % 2 == 1;   // nuclei 10LZZZAAAI: fermion for odd A
    return (id/1000) % 10 != 0;                     // baryons have a non-zero thousands digit
  }

  int MuCapBNVGenerator::baryonNumber(int pdgId) {
    const int id = std::abs(pdgId);
    if(id >= 1000000000) return ((pdgId > 0) ? 1 : -1) * ((id/10) % 1000);
    if(id >= 1000 && (id/1000) % 10 != 0) return (pdgId > 0) ? 1 : -1;
    return 0;
  }

  int MuCapBNVGenerator::leptonNumber(int pdgId) {
    const int id = std::abs(pdgId);
    if(id >= 11 && id <= 18) return (pdgId > 0) ? 1 : -1;
    return 0;
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::MuCapBNVGenerator)
