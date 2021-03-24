////////////////////////////////////////////////////////////////////////
// Class:       CRYGenPlots
// Plugin Type: analyzer (art v2_10_04)
// File:        CRYGenPlots_module.cc
//
// Generated at Wed Jun 27 18:05:40 2018 by Hoai Nam Tran using cetskelgen
// from cetlib version v3_02_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
namespace mu2e {
  class CRYGenPlots;
}

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;


class mu2e::CRYGenPlots : public art::EDAnalyzer {
  public:
    explicit CRYGenPlots(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRYGenPlots(CRYGenPlots const &) = delete;
    CRYGenPlots(CRYGenPlots &&) = delete;
    CRYGenPlots & operator = (CRYGenPlots const &) = delete;
    CRYGenPlots & operator = (CRYGenPlots &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:
    std::string processName_;
    std::string CRYModuleLabel_;
    std::string CRYInstanceName_;
    double _keMax;

    // histograms
    TH2F *_hXZ;
    TH1F *_hY;
    TH1F *_hKE;
    TH1F *_hTheta;
    TH1F *_hPhi;
    TH1F *_hPmag;
    TH1F *_hPyOverPmag;
    TH1F *_hTime;
    TH2F *_hPtypeKE;
    TH1F *_hNSecondaries;

    TTree *_cosmicTree;

    float _x;
    float _y;
    float _z;
    float _px;
    float _py;
    float _pz;
    float _theta;
    float _phi;
    float _KE;
    float _p;
    float _t;
    int _pdgId;

    void bookHists(art::ServiceHandle<art::TFileService> &);
    GlobalConstantsHandle<ParticleDataTable> pdt;
};


  mu2e::CRYGenPlots::CRYGenPlots(fhicl::ParameterSet const & p)
: EDAnalyzer(p)
  , processName_(p.get<std::string>("processName", ""))
  , CRYModuleLabel_(p.get<std::string>("CRYModuleLabel", "cryGen"))
  , CRYInstanceName_(p.get<std::string>("CRYInstanceName", ""))
  , _keMax(p.get<double>("keMax", 10E3))
{
  art::ServiceHandle<art::TFileService> tfs;

  bookHists(tfs);

  _cosmicTree = tfs->make<TTree>("cosmicTree", "TTree with cosmic ray info");

  _cosmicTree->Branch("x", &_x, "x/F");
  _cosmicTree->Branch("y", &_y, "y/F");
  _cosmicTree->Branch("z", &_z, "z/F");
  _cosmicTree->Branch("px", &_px, "px/F");
  _cosmicTree->Branch("py", &_py, "py/F");
  _cosmicTree->Branch("pz", &_pz, "pz/F");
  _cosmicTree->Branch("theta", &_theta, "theta/F");
  _cosmicTree->Branch("phi", &_phi, "phi/F");
  _cosmicTree->Branch("KE", &_KE, "KE/F");
  _cosmicTree->Branch("p", &_p, "p/F");
  _cosmicTree->Branch("t", &_t, "t/F");
  _cosmicTree->Branch("pdgId", &_pdgId, "pdgId/I");
}


void mu2e::CRYGenPlots::analyze(art::Event const & e)
{
  art::Handle<GenParticleCollection> gpHandle;
  bool success;

  if (processName_.length() > 0)
     success = e.getByLabel(CRYModuleLabel_, CRYInstanceName_, processName_,
         gpHandle);
  else
    success = e.getByLabel(CRYModuleLabel_, CRYInstanceName_, gpHandle);

  if (!success) 
    return;

  const auto & particles = *gpHandle;
  _hNSecondaries->Fill(particles.size());
  for(const auto & p : particles)
  {
    _hXZ->Fill(p.position().x(), p.position().z());
    _hY->Fill(p.position().y());

    const HepPDT::ParticleData& p_data = pdt->particle(p.pdgId()).ref();
    double mass = p_data.mass().value(); // in MeV

    HepLorentzVector mom4 = p.momentum();
    Hep3Vector mom3(mom4.px(), mom4.py(), mom4.pz());
    Hep3Vector yMom3(mom4.pz(), mom4.px(), -mom4.py()); // for theta, phi

    _hKE->Fill(mom4.e() - mass);
    _hTheta->Fill(yMom3.theta());
    _hPhi->Fill(yMom3.phi());

    _hPmag->Fill(mom3.mag());
    _hTime->Fill(p.time());
    _hPyOverPmag->Fill(mom4.py() / mom3.mag());

    _x = p.position().x();
    _y = p.position().y();
    _z = p.position().z();

    _px = mom4.px();
    _py = mom4.py();
    _pz = mom4.pz();

    _theta = yMom3.theta();
    _phi = yMom3.phi();

    _KE = mom4.e() - mass;
    _p = mom3.mag();

    _t = p.time();

    _pdgId = p.pdgId();
    _cosmicTree->Fill();

    _hKE->Fill(mom4.e() - mass);
    _hTheta->Fill(yMom3.theta());
    _hPhi->Fill(yMom3.phi());

    _hPmag->Fill(mom3.mag());
    _hTime->Fill(p.time());
    _hPyOverPmag->Fill(mom4.py() / mom3.mag());

    switch (p.pdgId()) {
      case 13: // mu-
        _hPtypeKE->Fill(mom4.e(), 0); break;
      case -13: // mu+
        _hPtypeKE->Fill(mom4.e(), 0); break;
      case 22: // photon
        _hPtypeKE->Fill(mom4.e(), 1); break;
      case -11: // e+
        _hPtypeKE->Fill(mom4.e(), 2); break;
      case 11: // e-
        _hPtypeKE->Fill(mom4.e(), 2); break;
      case 2112: // neutron
        _hPtypeKE->Fill(mom4.e(), 3); break;
      case -2112: // neutron
        _hPtypeKE->Fill(mom4.e(), 3); break;
      case 2212: // proton
        _hPtypeKE->Fill(mom4.e(), 4); break;
      case -2212: // proton
        _hPtypeKE->Fill(mom4.e(), 4); break;
      case 111: // pi0
        _hPtypeKE->Fill(mom4.e(), 5); break;
      case 211: // pi+
        _hPtypeKE->Fill(mom4.e(), 5); break;
      case -211: // pi-
        _hPtypeKE->Fill(mom4.e(), 5); break;
      case 130: // k0 L
        _hPtypeKE->Fill(mom4.e(), 6); break;
      case 310: // k0 S
        _hPtypeKE->Fill(mom4.e(), 6); break;
      case 311: // k0
        _hPtypeKE->Fill(mom4.e(), 6); break;
      case 321: // k+
        _hPtypeKE->Fill(mom4.e(), 6); break;
      case -321: // k-
        _hPtypeKE->Fill(mom4.e(), 6); break;
      default: // others
        _hPtypeKE->Fill(mom4.e(), 7); break;
    }

  }

}

void mu2e::CRYGenPlots::beginJob()
{
  // Implementation of optional member function here.
}

void mu2e::CRYGenPlots::endJob()
{
  // Implementation of optional member function here.
}

void mu2e::CRYGenPlots::bookHists(art::ServiceHandle<art::TFileService> &tfs)
{

  _hXZ = tfs->make<TH2F>("XZ", "XZ", 500, -2.0e5,  2.0e5, 500, -2.0e5, 2.0e5 );
  _hY = tfs->make<TH1F>("Y", "Y", 500, -15.0e3, 21.0e3 );

  _hKE = tfs->make<TH1F>("E", "E", 2000, 0, _keMax);
  _hTheta = tfs->make<TH1F>("Theta", "Theta",
      200, -M_PI - 0.5, M_PI + 0.5);
  _hPhi = tfs->make<TH1F>("Phi", "Phi",
      200, -M_PI - 0.5, M_PI + 0.5);

  _hPmag = tfs->make<TH1F>("Pmag", "Momentum modulus", 500, 0., 2E6);
  _hPyOverPmag = tfs->make<TH1F>("PyOverPmag", "Py/Pmag", 200, -20., 20.);
  _hTime = tfs->make<TH1F>("Time", "Timing of secondary particles",
      1000, -1.E-3, 1E0);

  _hNSecondaries = tfs->make<TH1F>("nSecondaries", "Number of secondaries",
      1000, 0, 1000);

  _hPtypeKE = tfs->make<TH2F>("PtypeKE", "Particle type vs energy",
      2000, 0., 2E6, 8, 0, 8);
  _hPtypeKE->GetYaxis()->SetBinLabel(1, "mu");
  _hPtypeKE->GetYaxis()->SetBinLabel(2, "gamma");
  _hPtypeKE->GetYaxis()->SetBinLabel(3, "electron");
  _hPtypeKE->GetYaxis()->SetBinLabel(4, "neutron");
  _hPtypeKE->GetYaxis()->SetBinLabel(5, "proton");
  _hPtypeKE->GetYaxis()->SetBinLabel(6, "pion");
  _hPtypeKE->GetYaxis()->SetBinLabel(7, "kaon");
  _hPtypeKE->GetYaxis()->SetBinLabel(8, "others");
}
DEFINE_ART_MODULE(mu2e::CRYGenPlots)
