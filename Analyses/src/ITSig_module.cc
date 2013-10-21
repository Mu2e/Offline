//
// A module to study background rates in the detector subsystems.
//
// $Id: ITSig_module.cc,v 1.3 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Gianni Onorato
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StatusG4.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include <cmath>
#include <iostream>
#include <set>
#include <memory>
#include <string>
#include <fstream>

using namespace std;

namespace mu2e {
  
  class ITSig : public art::EDAnalyzer {
  public:
    
    explicit ITSig(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _swiresStepPoints(pset.get<string>("swiresStepPoints","trackerSWires")),
      _fwiresStepPoints(pset.get<string>("fwiresStepPoints","itrackerFWires")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _extractElectronsData(pset.get<string>("elextractModuleLabel")),
      _minimumEnergyTracker(pset.get<double>("minimumEnergyTracker",0.0001)), // MeV
      _vdStepPoints(pset.get<string>("vdStepPoints","virtualdetector")),
      _nAnalyzed(0),
      _tNtup(0),
      _nBadG4Status(0),
      _nOverflow(0),
      _nKilled(0),
      _totalcputime(0),
      _totalrealtime(0)
    {
      cout << "Module ITSig is starting" << endl;
    }
    virtual ~ITSig() {
    }
    virtual void beginJob();
    virtual void endJob();
    
    void analyze(art::Event const& e );
    
    fstream* wirestxt;
    
  private:
    
    void doITracker(art::Event const& evt);
    
    // Diagnostic level
    int _diagLevel;
    
    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;
    std::string _swiresStepPoints;
    std::string _fwiresStepPoints;
    
    // Label of the module that made the hits.
    std::string _makerModuleLabel;
    
    // Label of the generator.
    std::string _generatorModuleLabel;
    
    // Label of the G4 module
    std::string _g4ModuleLabel;
    
    // Label of the module that made the hits.
    std::string _extractElectronsData;
    
    double _minimumEnergyTracker;  //minimum energy deposition of hits
    
    
    std::string _vdStepPoints;
    
    //number of analyzed events
    int _nAnalyzed;
    
    TNtuple* _tNtup;
    
    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;
    int _nVaribles;
    
  };
  
  void ITSig::beginJob( ) {
  }

  void ITSig::analyze(art::Event const& evt ) {
    
    ++_nAnalyzed;
    
    //*****test code******
    static int ncalls(0);
    ++ncalls;
    
    art::Handle<StatusG4> g4StatusHandle;
    evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    
    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
	<< "Aborting ITSig::analyze due to G4 status\n"
	<< g4Status;
      return;
    }
    
    if (g4Status.overflowSimParticles()) {
      ++_nOverflow;
      mf::LogError("G4")
	<< "Aborting ITSig::analyze due to overflow of particles\n"
	<< g4Status;
      return;
    }
    
    if (g4Status.nKilledStepLimit() > 0) {
      ++_nKilled;
      mf::LogError("G4")
	<< "Aborting ITSig::analyze due to nkilledStepLimit reached\n"
	<< g4Status;
      return;
    }
    
    _totalcputime += g4Status.cpuTime();
    _totalrealtime += g4Status.realTime();
    
    art::ServiceHandle<GeometryService> geom;
    
    if (ncalls == 1) {
      
      // cout << "This should be done only in the first event" << endl;
      
      wirestxt = new fstream("HitWires.txt",ios::out);
      
      art::ServiceHandle<art::TFileService> tfs;
      
      std::string generalInfo ="evt:run:nCells:nCellsSig:eDep:eDepSig:firstSigP:lastSigP:SigTotPath:SigInChamberTotPath";
      std::string electLoopInfo =":SigNloop:SigFrstLoopNHits:SigFrstLoopPath:SigScndLoopNHits:SigScndLoopPath:SigTrdLoopNHits:SigTrdLoopPath";
      std::string hitSWireInfo =":hitSWires:stepInSWires:pathInSWires:hitSWiresSig:stepInSWiresSig:pathInSWiresSig";
      std::string hitFWireInfo =":hitFWires:stepInFWires:pathInFWires:hitFWiresSig:stepInFWiresSig:pathInFWiresSig";
      std::string vdInfo =":enteringVd:enteringTime:enteringP:exitingVd:exitingTime:exitingP";
      _nVaribles = std::count(generalInfo.begin(), generalInfo.end(), ':')+1;
      _nVaribles += std::count(electLoopInfo.begin(), electLoopInfo.end(), ':');
      _nVaribles += std::count(hitSWireInfo.begin(), hitSWireInfo.end(), ':');
      _nVaribles += std::count(hitFWireInfo.begin(), hitFWireInfo.end(), ':');
      _nVaribles += std::count(vdInfo.begin(), vdInfo.end(), ':');
      _tNtup        = tfs->make<TNtuple>( "CellHits", "Cell Ntuple", (generalInfo+electLoopInfo+hitSWireInfo+hitFWireInfo+vdInfo).c_str() );
      
    }
    
    doITracker(evt);
    
  } // end of analyze
  
  void ITSig::endJob() {
    
    wirestxt->close();
    delete wirestxt;
    
    cout << "ITSig::endJob Number of events skipped "
	 << "due to G4 completion status: "
	 << _nBadG4Status
	 << "\nITSig::endJob Number of overflow events "
	 << "due to too many particles in G4: "
	 << _nOverflow
	 << "\nITSig::endJob Number of events with killed particles "
	 << "due to too many steps in G4: "
	 << _nKilled
	 << "\nITSig::endJob total CpuTime "
	 << _totalcputime
	 << "\nITSig::endJob total RealTime "
	 << _totalrealtime
	 << endl;
  }
  
  void ITSig::doITracker(art::Event const& evt) {
    
    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();
    
    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
    
    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);
    
    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);
    
    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);
    

    //Handle to VD steps
    art::Handle<StepPointMCCollection> vdHits;
    evt.getByLabel(_g4ModuleLabel,_vdStepPoints,vdHits);

    // Find original G4 steps in the wires;
    art::Handle<StepPointMCCollection> wShits;
    evt.getByLabel(_g4ModuleLabel,_swiresStepPoints,wShits);
    art::Handle<StepPointMCCollection> wFhits;
    evt.getByLabel(_g4ModuleLabel,_fwiresStepPoints,wFhits);
    
    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    evt.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
    
    // cout << "event " << evt.id().event() << ": reading stuff" << endl;

    size_t nSWHits = wShits->size();
    size_t nFWHits = wFhits->size();
    
    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );
    
    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }
    
    size_t nCells = hits->size();
    
    float nCellsSig(0);
    double eDep(0), eDepSig(0);
    
    float *tntpArray = new float [_nVaribles];
    for (int iv=0; iv<_nVaribles; ++iv) {
      tntpArray[iv]=0.0;
    }

    // cout << "event " << evt.id().event() << ": first fills" << endl;

    int idx(0);
    tntpArray[idx++] = evt.id().event(); //1
    tntpArray[idx++] = evt.run();//2
    
    SimParticleCollection::key_type trackId;
    // cout << "event " << evt.id().event() << ": looking at cells" << endl;
    
    double tempsigfirsttime = 10e9;
    int jsigfirst = -1;
    int isigfirst = -1;
    double tempsiglasttime = -100;
    int jsiglast = -1;
    int isiglast = -1;

    for (size_t i=0; i<nCells; ++i) {
      // Access data
      StrawHit             const&  hit(hits->at(i));
      PtrStepPointMCVector const&  mcptr(hits_mcptr->at(i));
      
      //Skip the straw if the energy of the hit is smaller than the minimum required
      if (hit.energyDep() < _minimumEnergyTracker) continue;
      
      eDep += hit.energyDep();
      
      //Find the first stepPointMC by time of flight
      double temptime = 10e9;
      size_t jfirst = 0;
      for (size_t j = 0; j < mcptr.size(); ++j) {
	if (mcptr[j]->time() < temptime) {
	  temptime = mcptr[j]->time();
	  jfirst = j;
	}
      }

      //Find the first sigstepPointMC by time of flight
      for (size_t j = 0; j < mcptr.size(); ++j) {
	if (simParticles->at(mcptr[j]->trackId()).fromGenerator()) {
	  if (mcptr[j]->time() < tempsigfirsttime) {
	    tempsigfirsttime = mcptr[j]->time();
	    jsigfirst = j;
	    isigfirst = i;
	  }
	}
      }

      //Find the last stepPointMC by time of flight
      for (size_t j = 0; j < mcptr.size(); ++j) {
	if (simParticles->at(mcptr[j]->trackId()).fromGenerator()) {
	  if (mcptr[j]->time() > tempsiglasttime) {
	  tempsiglasttime = mcptr[j]->time();
	  jsiglast = j;
	  isiglast = i;
	  }
	}
      }

      // cout << "cell " << i << ": looking at first steppoint" << endl;

      StepPointMC const& mchit = *mcptr[jfirst];
      
      // The simulated particle that made this hit.
      trackId=mchit.trackId();
      SimParticle const& sim = simParticles->at(trackId);
      if ( sim.fromGenerator() && sim.pdgId()==11 ) {
	eDepSig += hit.energyDep();
	nCellsSig += 1;
      }
    }
    
    tntpArray[idx++] = nCells;
    tntpArray[idx++] = nCellsSig;
    tntpArray[idx++] = eDep;
    tntpArray[idx++] = eDepSig;
    
    if (jsigfirst != -1) {
      PtrStepPointMCVector const&  mcptrFirstV(hits_mcptr->at(isigfirst));
      PtrStepPointMCVector const&  mcptrLastV(hits_mcptr->at(isiglast));
      tntpArray[idx++] = mcptrFirstV[jsigfirst]->momentum().mag();
      tntpArray[idx++] = mcptrLastV[jsiglast]->momentum().mag();
    } else {
      tntpArray[idx++] = -10;
      tntpArray[idx++] = -10;
    }
    if (genEltrks->size()>0) {

      for ( std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
	VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
	if (iEltrk.getTrkID()==trackId) {
	  
	  //	  cout << "element 0 done. Go with the 1" << endl;
	  tntpArray[idx++] = iEltrk.getTotPath();
	  //	  cout << "element 1 done. Go with the 2" << endl;
	  tntpArray[idx++] = iEltrk.getInChamberTotPath();
	  //	  cout << "element 2 done. Go with the 3" << endl;
	  tntpArray[idx++] = iEltrk.getNumOfLoops();
	  // cout << "element 3 done. Go with the 4" << endl;
	  // cout << "element 4 done. Go with the 5" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopNHit(0);
	  // cout << "element 5 done. Go with the 6" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopTOF(0)*CLHEP::c_light;
	  // cout << "element 6 done. Go with the 7" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopNHit(1);
	  // cout << "element 7 done. Go with the 8" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopTOF(1)*CLHEP::c_light;
	  // cout << "element 8 done. Go with the 9" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopNHit(2);
	  // cout << "element 9 done. Go with the 10" << endl;
	  tntpArray[idx++] = iEltrk.getithLoopTOF(2)*CLHEP::c_light;
	  // cout << "Cool!" << endl;
	}
      }
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
    }

    //    float nHitWires(0);

    // cout << "event " << evt.id().event() << ": looking at wires" << endl;    
    if (nSWHits>0) {
      
      //      cout << "Some sense wire has been hit!" << endl;
      
      size_t wSSteps = wShits->size();
      
      //      cout << "Precisely they are " << wSSteps << endl;
      
      float wSStepsSig(0);
      
      set<unsigned long> wires, sigwires;
      double totPathInWires(0.0), totPathInWiresSig(0.0);
    
      for (size_t i=0; i<wSSteps; ++i) {
	// cout << "step " << i << ": looking at sense wire steps" << endl;
	StepPointMC const& whit = (*wShits)[i];
	SimParticleCollection::key_type trackId = whit.trackId();
	SimParticle const* sim = simParticles->getOrNull(trackId);
	if( !sim ) continue;
	wires.insert(whit.volumeId());
	totPathInWires+=whit.stepLength();
	*wirestxt << "Hit n. " << i+1 << " in the sense wire " << whit.volumeId()
		  << " made by " << sim->pdgId() << endl;
	if (sim->fromGenerator()) {
	  wSStepsSig+=1;
	  sigwires.insert(whit.volumeId());
	  totPathInWiresSig+=whit.stepLength();
	}
      }
      tntpArray[idx++] = wires.size();
      tntpArray[idx++] = nSWHits;
      tntpArray[idx++] = totPathInWires;
      tntpArray[idx++] = sigwires.size();
      tntpArray[idx++] = wSStepsSig;
      tntpArray[idx++] = totPathInWiresSig;
      //_tNtup->Fill(tntpArray);
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
    }
    
    
    if (nFWHits>0) {
      
      //      cout << "Some field wire has been hit!" << endl;
      
      size_t wFSteps = wFhits->size();
      
      //      cout << "Precisely they are " << wFSteps << endl;
      
      float wFStepsSig(0);
      
      set<unsigned long> wires, sigwires;
      double totPathInWires(0.0), totPathInWiresSig(0.0);
      
      for (size_t i=0; i<wFSteps; ++i) {
	// cout << "step " << i << ": looking at field wire steps" << endl;
	StepPointMC const& whit = (*wFhits)[i];
	SimParticleCollection::key_type trackId = whit.trackId();
	SimParticle const* sim = simParticles->getOrNull(trackId);
	if( !sim ) continue;
	wires.insert(whit.volumeId());
	totPathInWires+=whit.stepLength();
	*wirestxt << "Hit n. " << i+1 << " in the field wire " << whit.volumeId()
		  << " made by " << sim->pdgId() << endl;
	if (sim->fromGenerator()) {
	  wFStepsSig+=1;
	  sigwires.insert(whit.volumeId());
	  totPathInWiresSig+=whit.stepLength();
	}
      }
      tntpArray[idx++] = wires.size();
      tntpArray[idx++] = nFWHits;
      tntpArray[idx++] = totPathInWires;
      tntpArray[idx++] = sigwires.size();
      tntpArray[idx++] = wFStepsSig;
      tntpArray[idx++] = totPathInWiresSig;
      //_tNtup->Fill(tntpArray);
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
    }
    
    
    //Loop over vd hits
    bool foundin = false;
    bool foundout = false;
    int vdIn(0), vdOut(0);
    double pIn(0), pOut(0);
    double intime = 1e15; 
    double outtime = -1e15;
    size_t inIndex(0), outIndex(0);
    // cout << "event " << evt.id().event() << ": looking at VD steps" << endl;
    if (vdHits.isValid()) {
      for (size_t i=0; i<vdHits->size(); ++i) {
	const StepPointMC& hit = (*vdHits)[i];
	if (simParticles->at(hit.trackId()).fromGenerator()) {
	  int id = hit.volumeId();
	  //if (id >= 40 && id <= 42) {
	  if (id == VirtualDetectorId::IT_VD_InSurf||
              id==VirtualDetectorId::IT_VD_EndCap_Front||
              id==VirtualDetectorId::IT_VD_EndCap_Back) {
	    if (hit.time() < intime) {
	      pIn = hit.momentum().mag();
	      foundin = true;
	      vdIn = id;
	      inIndex = i;
	      intime = hit.time();
	    }
	    if (hit.time() > outtime) {
	      pOut = hit.momentum().mag();
	      foundout = true;
	      vdOut = id;
	      outIndex = i;
	      outtime = hit.time();
	    }
	  }
	}
      }
    }
   
    if (foundin && foundout && (inIndex!=outIndex) && (intime < outtime)) {
      tntpArray[idx++] = vdIn;
      tntpArray[idx++] = intime;
      tntpArray[idx++] = pIn;
      tntpArray[idx++] = vdOut;
      tntpArray[idx++] = outtime;
      tntpArray[idx++] = pOut;
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
    }
    
    
    
    
    _tNtup->Fill(tntpArray);
    
    delete [] tntpArray;
    
  } // end of doITracker
  
}

using mu2e::ITSig;
DEFINE_ART_MODULE(ITSig);

