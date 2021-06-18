////////////////////////////////////////////////////////////////////////
// Class:       CORSIKAGenPlots
// Plugin Type: analyzer (art v2_10_04)
// File:        CORSIKAGenPlots_module.cc
//
// Stefano Roberto Soleti (roberto@lbl.gov)
//
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
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
namespace mu2e {
  class CORSIKAGenPlots;
}

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;


class mu2e::CORSIKAGenPlots : public art::EDAnalyzer {
  public:
    explicit CORSIKAGenPlots(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CORSIKAGenPlots(CORSIKAGenPlots const &) = delete;
    CORSIKAGenPlots(CORSIKAGenPlots &&) = delete;
    CORSIKAGenPlots & operator = (CORSIKAGenPlots const &) = delete;
    CORSIKAGenPlots & operator = (CORSIKAGenPlots &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:
    std::string processName_;
    std::string CORSIKAModuleLabel_;
    std::string CORSIKAInstanceName_;
    float _keMax = std::numeric_limits<float>::max();

    // histograms
    TH2F *_hXZ = nullptr;
    TH1F *_hY = nullptr;
    TH1F *_hKE = nullptr;
    TH1F *_hTheta = nullptr;
    TH1F *_hPhi = nullptr;
    TH1F *_hPmag = nullptr;
    TH1F *_hPyOverPmag = nullptr;
    TH1F *_hTime = nullptr;
    TH2F *_hPtypeKE = nullptr;
    TH1F *_hNSecondaries = nullptr;

    TTree *_cosmicTree = nullptr;
    TTree *_eventTree = nullptr;

    float _x = std::numeric_limits<float>::lowest();
    float _y = std::numeric_limits<float>::lowest();
    float _z = std::numeric_limits<float>::lowest();
    float _px = std::numeric_limits<float>::lowest();
    float _py = std::numeric_limits<float>::lowest();
    float _pz = std::numeric_limits<float>::lowest();
    float _theta = std::numeric_limits<float>::lowest();
    float _phi = std::numeric_limits<float>::lowest();
    float _KE = std::numeric_limits<float>::lowest();
    float _p = std::numeric_limits<float>::lowest();
    float _t = std::numeric_limits<float>::lowest();
    std::vector<float> _dca;
    PDGCode::type _pdgId;

    void bookHists(art::ServiceHandle<art::TFileService> &);
    GlobalConstantsHandle<ParticleDataTable> pdt;
};

mu2e::CORSIKAGenPlots::CORSIKAGenPlots(fhicl::ParameterSet const &p)
    : EDAnalyzer(p),
    processName_(p.get<std::string>("processName", "")),
    CORSIKAModuleLabel_(p.get<std::string>("CORSIKAModuleLabel", "FromCorsikaBinary")),
    CORSIKAInstanceName_(p.get<std::string>("CORSIKAInstanceName", "")),
    _keMax(p.get<double>("keMax", 10E3))
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

  _eventTree = tfs->make<TTree>("eventTree", "TTree with cosmic ray info per shower");
  _eventTree->Branch("dca", &_dca);
}


void mu2e::CORSIKAGenPlots::analyze(art::Event const & e)
{
  art::Handle<GenParticleCollection> gpHandle;
  bool success;

  if (processName_.length() > 0)
     success = e.getByLabel(CORSIKAModuleLabel_, CORSIKAInstanceName_, processName_,
         gpHandle);
  else
    success = e.getByLabel(CORSIKAModuleLabel_, CORSIKAInstanceName_, gpHandle);

  if (!success)
    return;

  const auto & particles = *gpHandle;

  // Store the point of closest approach between target box and roof for
  // events with more than one particle
  _dca.clear();
  for (GenParticleCollection::const_iterator i = particles.begin(); i != particles.end(); ++i)
  {
    for (GenParticleCollection::const_iterator j = i+1; j != particles.end(); ++j)
    {
        GenParticle const &particle1 = *i;
        GenParticle const &particle2 = *j;

        TwoLinePCA twoLine(particle1.position(), -particle1.momentum().vect(), particle2.position(), -particle2.momentum().vect());
        const CLHEP::Hep3Vector point1 = twoLine.point1();
        const CLHEP::Hep3Vector point2 = twoLine.point2();


        if (point1.y() > 5000 &&
            point1.y() < 15365.4 &&
            point2.y() > 5000 &&
            point2.y() < 15365.4)
        {
          _dca.push_back(twoLine.dca());
        }
    }
  }
  _eventTree->Fill();

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

void mu2e::CORSIKAGenPlots::beginJob()
{
  // Implementation of optional member function here.
}

void mu2e::CORSIKAGenPlots::endJob()
{
  // Implementation of optional member function here.
}

void mu2e::CORSIKAGenPlots::bookHists(art::ServiceHandle<art::TFileService> &tfs)
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
DEFINE_ART_MODULE(mu2e::CORSIKAGenPlots)
