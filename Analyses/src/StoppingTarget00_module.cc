//
// A first look at muons stopping in stopping targets.
//
// $Id: StoppingTarget00_module.cc,v 1.14 2013/10/21 21:15:46 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:15:46 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"

// Mu2e includes.
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"

// Root includes.
#include "TH1F.h"
#include "TNtuple.h"

using namespace std;

namespace mu2e {

  class StoppingTarget00 : public art::EDAnalyzer {
  public:
    explicit StoppingTarget00(fhicl::ParameterSet const& pset);
    virtual ~StoppingTarget00() { }

    void beginJob();
    void endJob();

    void beginRun(art::Run const& );
    void analyze( art::Event const& );

  private:

    // Label of the module that created the data products.
    std::string _g4ModuleLabel;

    // Instance names of data products
    std::string _targetStepPoints;
    std::string _vdStepPoints;

    // Name of file to hold get the points and times that muons stop.
    std::string _muonPointFile;

    // Offset to put coordinates in a special reference frame:
    //  - (x,y) origin on DS axis; z origin at Mu2e origin.
    CLHEP::Hep3Vector _dsOffset;

    TH1F* _hStopFoil;
    TH1F* _hnSimPart;
    TH1F* _hnSimPartZoom;
    TH1F* _hzInFoil;
    TNtuple* _nt;

    map<int,int> _nonMuon;
    map<int,int> _particleCounts;
    map<int,int> _vdGlobalCounts;

    map<ProcessCode,int> _stopCodes;
    map<PhysicalVolumeInfo const*,int> _startVols;
    map<PhysicalVolumeInfo const*,int> _stoppingVols;

  };

  StoppingTarget00::StoppingTarget00(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel"))
    ,_targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget"))
    ,_vdStepPoints(pset.get<string>("vdStepPoints","virtualdetector"))
    ,_muonPointFile(pset.get<string>("muonPointFile",""))
    ,_dsOffset()
    ,_hStopFoil(0)
    ,_hnSimPart(0)
    ,_hnSimPartZoom(0)
    ,_hzInFoil(0)
    ,_nt(0)
    ,_nonMuon()
    ,_particleCounts()
    ,_stopCodes()
  {}

  void StoppingTarget00::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    _hStopFoil = tfs->make<TH1F>( "hStopFoil", "Number of Stopping Foil;(mm)", 34, 0., 17. );
    _hnSimPart = tfs->make<TH1F>( "hnSimPart", "Number of SimParticles", 200, 0., 1000. );
    _hnSimPartZoom = tfs->make<TH1F>( "hnSimPartZoom", "Number of SimParticles", 50, 0., 50. );
    _hzInFoil  = tfs->make<TH1F>( "hzInFoil", "z Position Within Foil(units of halfThickness", 
                                  100, -1., 1. );
    _nt = tfs->make<TNtuple>( "nt", "Muon ntuple",
                                "cx:cy:cz:cp:cpt:cpz:cke:ct:sx:sy:sz:sp:spt:st:stau:scode:edep:sfoil:nfoils:x9:y9:pt9:pz9:x10:y10:pt10:pz10");

  }

  void StoppingTarget00::beginRun(art::Run const& run ){

    // Information about the detector coordinate system.
    GeomHandle<DetectorSystem> det;
    _dsOffset = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );

    /*
    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    run.getByLabel(_g4ModuleLabel, volsHandle);
    PhysicalVolumeInfoCollection const& vols(*volsHandle);


    for ( size_t i = 0; i<vols.size(); ++i){
      cout << "Info: "
           << i  << " "
           << vols.at(i)
           << endl;
    }
    */
  }

  void
  StoppingTarget00::analyze(art::Event const& event ) {

    // Information about the detector coordinate system.
    //GeomHandle<DetectorSystem> det;

    GeomHandle<StoppingTarget> target;

    // Simulated particles.
    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(_g4ModuleLabel,simsHandle);
    SimParticleCollection const& sims(*simsHandle);
    if ( sims.size() == 0 ){
      mf::LogInfo("G4")
        << "No particles in SimParticleCollection.  Hope that's OK.";
      return;
    }

    // Steps in the stopping target
    art::Handle<StepPointMCCollection> stHitsHandle;
    event.getByLabel(_g4ModuleLabel,_targetStepPoints,stHitsHandle);
    StepPointMCCollection const& sthits(*stHitsHandle);

    // Steps in the virtual detectors.
    art::Handle<StepPointMCCollection> vdhitsHandle;
    event.getByLabel(_g4ModuleLabel,_vdStepPoints,vdhitsHandle);
    StepPointMCCollection const& vdhits(*vdhitsHandle);

    // Information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    event.getRun().getByLabel(_g4ModuleLabel, volsHandle);
    PhysicalVolumeInfoCollection const& vols(*volsHandle);

    _hnSimPart->Fill(sims.size());
    _hnSimPartZoom->Fill(sims.size());

    // In my files the muon is always simparticle 1.
    SimParticleCollection::key_type muIndex(1);
    SimParticle const& mu = sims.at(muIndex);
    if ( mu.pdgId() != PDGCode::mu_minus ) {
      ++_nonMuon[mu.pdgId()];
      return;
    }

    PhysicalVolumeInfo const& startVol = vols.at(mu.startVolumeIndex());
    PhysicalVolumeInfo const&   endVol = vols.at(mu.endVolumeIndex());

    ++_stopCodes[mu.stoppingCode()];
    ++_stoppingVols[&endVol];
    ++_startVols[&startVol];

    int stopFoilId(-1);
    if ( endVol.name().compare(0,11,"TargetFoil_") == 0 ) {
      stopFoilId = endVol.copyNo();
      _hStopFoil->Fill(stopFoilId);
    }

    // Sum of energy loss in the foils.
    double sumEDep(0.);

    set<int> hitFoils;
    int nFoils(0);
    for ( StepPointMCCollection:: const_iterator i=sthits.begin();
          i !=sthits.end(); ++i ){
      StepPointMC const& hit = *i;
      if ( hit.trackId() != muIndex ) continue;
      sumEDep += hit.eDep();
      hitFoils.insert(hit.volumeId());
      ++nFoils;
    }

    /*
    if ( !hitFoils.empty() && nFoils > (int(hitFoils.size())+1) ){
      cout << "Mark: "
           << event.id() << " "
           << endVol << " "
           << nFoils << " "
           << hitFoils.size() << " "
           << muIndex
           << endl;
      int nn(0);
      set<int> dummy;
      for ( StepPointMCCollection:: const_iterator i=sthits.begin();
            i !=sthits.end(); ++i ){
        StepPointMC const& hit = *i;
        if ( hit.trackId() != muIndex ) continue;
        dummy.insert(hit.volumeId());
        cout << "Hit Foil: "
             << nn++ <<  " "
             << hit.trackId() << " "
             << hit.volumeId() << " | "
             << dummy.size() << " | "
             << hit.time() << " "
             << hit.properTime() <<  " "
             << hit.stepLength() <<  " "
             << hit.totalEDep() <<  " | "
             << hit.momentum() <<  " | "
             << hit.momentum().mag()
             << endl;
      }
    }
    */

    struct ntInfo{
      float cx;  // Position, momentum, kinetic energy and time at start
      float cy;
      float cz;
      float cp;
      float cpt;
      float cpz;
      float cke;
      float ct;
      float ex;  // Same at end
      float ey;
      float ez;
      float ep;
      float ept;
      float et;
      float tau;   // Proper time at end
      float ecode; // Reason why it stoppd
      float eDep;
      float sFoil;
      float nFoils;
      float x9;     // Position at virtual detector 9
      float y9;
      float pt9;
      float pz9;
      float x10;    // Position at virtual detector 10
      float y10;
      float pt10;
      float pz10;
      ntInfo():
        cx(0), cy(0), cz(0), cp(0), cke(0), ct(0),
        ex(0), ey(0), ez(0), ep(0),         et(0),
        tau(0),
        ecode(0),
        eDep(0),
        nFoils(0),
        x9(0),   y9(0),
        x10(0), y10(0){}
    }xx;

    set<int> vdCounts;
    int nvd(0);
    for ( StepPointMCCollection:: const_iterator i=vdhits.begin();
          i !=vdhits.end(); ++i ){
      StepPointMC const& hit = *i;
      if ( hit.trackId() != muIndex ) continue;

      // Hits are in the Mu2e coordinate system.  Report these in detector system.
      if ( hit.volumeId() == 9 ){
        CLHEP::Hep3Vector pos = hit.position() - _dsOffset;
        xx.x9  = pos.x();
        xx.y9  = pos.y();
        xx.pt9 = hit.momentum().perp();
        xx.pz9 = hit.momentum().z();
      }
      if ( hit.volumeId() == 10 ){
        CLHEP::Hep3Vector pos = hit.position() - _dsOffset;
        xx.x10  = pos.x();
        xx.y10  = pos.y();
        xx.pt10 = hit.momentum().perp();
        xx.pz10 = hit.momentum().z();
      }
      /*
      cout << "VD: "
           << hit.volumeId() << " "
           << hit.position().z() << " "
           << hit.momentum().mag() << " "
           << hit.eDep()*1000.
           << endl;
      */
      vdCounts.insert(hit.volumeId());
      ++_vdGlobalCounts[hit.volumeId()];
      ++nvd;
    }

    /*
    double deTest = mu.startMomentum().e()-mu.endMomentum().e() - sumEDep;
    cout << "Energetics: "
         << mu.startMomentum().e() << " "
         << mu.endMomentum().e() << " "
         << mu.startMomentum().e()-mu.endMomentum().e() << " "
         << sumEDep << " "
         << deTest << " "
         << mu.stoppingCode() << " "
         << endVol << " | "
         << nFoils << " "
         << hitFoils.size() << " | "
         << vdCounts.size() << " "
         << nvd
         << endl;
    */

    // Transform point from Mu2e frame to detector frame.
    CLHEP::Hep3Vector startPos = mu.startPosition() - _dsOffset;
    CLHEP::Hep3Vector endPos   = mu.endPosition()   - _dsOffset;

    xx.cx     = startPos.x();
    xx.cy     = startPos.y();
    xx.cz     = startPos.z();
    xx.cp     = mu.startMomentum().vect().mag();
    xx.cpt    = mu.startMomentum().vect().perp();
    xx.cpz    = mu.startMomentum().vect().z();
    xx.cke    = mu.startMomentum().e()-mu.startMomentum().mag();
    xx.ct     = mu.startGlobalTime();
    xx.ex     = endPos.x();
    xx.ey     = endPos.y();
    xx.ez     = endPos.z();
    xx.ep     = mu.endMomentum().vect().mag();
    xx.ept    = mu.endMomentum().vect().perp();
    xx.et     = mu.endGlobalTime();
    xx.tau    = mu.endProperTime();
    xx.ecode  = mu.stoppingCode();
    xx.eDep   = sumEDep;
    xx.sFoil  = stopFoilId;
    xx.nFoils = hitFoils.size();

    _nt->Fill( &xx.cx);

    if ( stopFoilId > -1 && mu.stoppingCode() == ProcessCode::muMinusCaptureAtRest ){
      TargetFoil const& foil = target->foil(stopFoilId);

      double dz = endPos.z() - foil.centerInDetectorSystem().z();
      _hzInFoil->Fill( dz/foil.halfThickness() );
      /*
      cout << "Test: "
           << endPos.z() << " "
           << foil.center().z() <<  " "
           << dz << " " 
           << dz/foil.halfThickness()
           << endl;
      */
      if ( !_muonPointFile.empty() ) {
        static ofstream fout(_muonPointFile.c_str());
        static bool first(true);
        if ( first ){
          fout << "begin data" << endl;
          first = false;
        }
        // Print out in Mu2e coordinates
        fout  << setprecision(8)
              << endPos.x() + _dsOffset.x() << " "
              << endPos.y() + _dsOffset.y() << " "
              << endPos.z() + _dsOffset.z() << " "
              << mu.endGlobalTime()
              << endl;
      }
    }
                               

    // Table of reasons why particles stopped.
    for ( SimParticleCollection::const_iterator i=sims.begin();
          i!=sims.end(); ++i ){
      SimParticle const& sim = i->second;
      ++_particleCounts[sim.pdgId()];
    }

  } // end of ::analyze.

  void StoppingTarget00::endJob(){
    cout << "\nNumber of distinct tracked particles: " << _particleCounts.size() << endl;
    int n(0);
    for ( map<int,int>::const_iterator i=_particleCounts.begin();
          i != _particleCounts.end(); ++i ){
      cout << "   Tracked Particle: "
           << n++ << " "
           << i->first << ": "
           << i->second
           << endl;
    }

    cout << "\nNumber of times first particle was not a muon: " << _nonMuon.size() << endl;
    n=0;
    for ( map<int,int>::const_iterator i=_nonMuon.begin();
          i != _nonMuon.end(); ++i ){
      cout << "   nonMuon: "
           << n++ << " "
           << i->first << ": "
           << i->second
           << endl;
    }

    cout << "\nNumber of unique stopping Codes: " << _stopCodes.size() << endl;
    n=0;
    for ( map<ProcessCode,int>::const_iterator i=_stopCodes.begin();
          i != _stopCodes.end(); ++i ){
      cout << "   Stopping Code: "
           << n++ << " "
           << i->first << ": "
           << i->second
           << endl;
    }

    cout << "\nNumber of different stopping volumes: " << _stoppingVols.size() << endl;
    n=0;
    for ( map<PhysicalVolumeInfo const*,int>::const_iterator i=_stoppingVols.begin();
          i != _stoppingVols.end(); ++i ){
      cout << "   Stopping Volume: "
           << n++ << " "
           << *i->first << ": "
           << i->second
           << endl;
    }

    cout << "\nNumber of different start volumes: " << _startVols.size() << endl;
    n=0;
    for ( map<PhysicalVolumeInfo const*,int>::const_iterator i=_startVols.begin();
          i != _startVols.end(); ++i ){
      cout << "   Start Volume: "
           << n++ << " "
           << *i->first << ": "
           << i->second
           << endl;
    }

    cout << "\nNumber of different virtual detectors: " << _vdGlobalCounts.size() << endl;
    n=0;
    for ( map<int,int>::const_iterator i=_vdGlobalCounts.begin();
          i != _vdGlobalCounts.end(); ++i ){
      cout << "   Virtual Detector: "
           << n++ << " "
           << i->first << ": "
           << i->second
           << endl;
    }


  }
}

using mu2e::StoppingTarget00;
DEFINE_ART_MODULE(StoppingTarget00);
