//
// A module to study background rates in the detector subsystems.
//
// $Id: TSig_module.cc,v 1.3 2014/09/03 15:51:18 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/03 15:51:18 $
//
// Original author G. Tassielli
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
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
#include <iterator>
#include <utility>

using namespace std;

namespace mu2e {
  
  class TSig : public art::EDAnalyzer {
  public:
    
    explicit TSig(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _swiresStepPoints(pset.get<string>("swiresStepPoints","trackerSWires")),
      _wallsStepPoints(pset.get<string>("wallsStepPoints","trackerWalls")),
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
      cout << "Module TSig is starting" << endl;
    }
    virtual ~TSig() {
    }
    virtual void beginJob();
    virtual void endJob();
    
    void analyze(art::Event const& e );
    
    fstream* wirestxt;
    
  private:
    
    void doTracker(art::Event const& evt);
    
    // Diagnostic level
    int _diagLevel;
    
    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;
    std::string _swiresStepPoints;
    std::string _wallsStepPoints;
    
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
    TTree*   _tSingleWallData;
    int      _nHit;
    double   _wallEloss[1000];
    double   _wallPath[1000];
    TTree*   _tSingleCellData;
    int      _nCells;
    double   _cellEloss[1000];
    double   _cellPath[1000];

    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;
    int _nVaribles;
    
  };
  
  void TSig::beginJob( ) {
  }

  void TSig::analyze(art::Event const& evt ) {
    
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
	<< "Aborting TSig::analyze due to G4 status\n"
	<< g4Status;
      return;
    }
    
    if (g4Status.overflowSimParticles()) {
      ++_nOverflow;
      mf::LogError("G4")
	<< "Aborting TSig::analyze due to overflow of particles\n"
	<< g4Status;
      return;
    }
    
    if (g4Status.nKilledStepLimit() > 0) {
      ++_nKilled;
      mf::LogError("G4")
	<< "Aborting TSig::analyze due to nkilledStepLimit reached\n"
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
      std::string hitSWireInfo =":hitSWires:stepInSWires:pathInSWires:eDepInSWires:hitSWiresSig:stepInSWiresSig:pathInSWiresSig:eDepInSWiresSig";
      std::string hitWallInfo =":hitWalls:stepInWalls:pathInWalls:eDepInWalls:hitWallsSig:stepInWallsSig:pathInWallsSig:eDepInWallsSig";
      std::string vdInfo =":enteringVd:enteringTime:enteringP:exitingVd:exitingTime:exitingP";
      _nVaribles = std::count(generalInfo.begin(), generalInfo.end(), ':')+1;
      _nVaribles += std::count(electLoopInfo.begin(), electLoopInfo.end(), ':');
      _nVaribles += std::count(hitSWireInfo.begin(), hitSWireInfo.end(), ':');
      _nVaribles += std::count(hitWallInfo.begin(), hitWallInfo.end(), ':');
      _nVaribles += std::count(vdInfo.begin(), vdInfo.end(), ':');
      _tNtup        = tfs->make<TNtuple>( "CellHits", "Cell Ntuple", (generalInfo+electLoopInfo+hitSWireInfo+hitWallInfo+vdInfo).c_str() );

      _tSingleWallData = tfs->make<TTree>( "singleWallData", "Eloss for each Cell");
      _tSingleWallData->Branch("nHit",&_nHit,"_nHit/I");
      _tSingleWallData->Branch("wallEloss",_wallEloss,"_wallEloss[_nHit]/D");
      _tSingleWallData->Branch("wallPath",_wallPath,"_wallPath[_nHit]/D");

      _tSingleCellData = tfs->make<TTree>( "singleCellData", "Eloss for each Cell");
      _tSingleCellData->Branch("nCells",&_nCells,"_nCells/I");
      _tSingleCellData->Branch("cellEloss",_cellEloss,"_cellEloss[_nCells]/D");
      _tSingleCellData->Branch("cellPath",_cellPath,"_cellPath[_nCells]/D");

    }
    
    doTracker(evt);
    
  } // end of analyze
  
  void TSig::endJob() {
    
    wirestxt->close();
    delete wirestxt;
    
    cout << "TSig::endJob Number of events skipped "
	 << "due to G4 completion status: "
	 << _nBadG4Status
	 << "\nTSig::endJob Number of overflow events "
	 << "due to too many particles in G4: "
	 << _nOverflow
	 << "\nTSig::endJob Number of events with killed particles "
	 << "due to too many steps in G4: "
	 << _nKilled
	 << "\nTSig::endJob total CpuTime "
	 << _totalcputime
	 << "\nTSig::endJob total RealTime "
	 << _totalrealtime
	 << endl;
  }
  
  void TSig::doTracker(art::Event const& evt) {
    
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
    art::Handle<StepPointMCCollection> wllhits;
    evt.getByLabel(_g4ModuleLabel,_wallsStepPoints,wllhits);
    
    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    evt.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
    
    // cout << "event " << evt.id().event() << ": reading stuff" << endl;

    size_t nSWHits = wShits->size();
    size_t nWlHits = wllhits->size();
    
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

    bool *sigCellNotDone = new bool [nCells];
    for (size_t iscd=0; iscd<nCells; ++iscd) { sigCellNotDone[iscd]=true; }

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

      if (sigCellNotDone[(int)nCellsSig]) {
              _cellEloss[(int)nCellsSig]=0;
              _cellPath[(int)nCellsSig]=0;
      }
      //Find the first sigstepPointMC by time of flight
      for (size_t j = 0; j < mcptr.size(); ++j) {
	if (simParticles->at(mcptr[j]->trackId()).fromGenerator() && simParticles->at(mcptr[j]->trackId()).pdgId()==11) {
	  if (mcptr[j]->time() < tempsigfirsttime) {
	    tempsigfirsttime = mcptr[j]->time();
	    jsigfirst = j;
	    isigfirst = i;
	  }

	  _cellEloss[(int)nCellsSig]+=mcptr[j]->eDep();
	  _cellPath[(int)nCellsSig]+=mcptr[j]->stepLength();
	  sigCellNotDone[(int)nCellsSig] = false;

	}
      }

      //Find the last stepPointMC by time of flight
      for (size_t j = 0; j < mcptr.size(); ++j) {
	if (simParticles->at(mcptr[j]->trackId()).fromGenerator() && simParticles->at(mcptr[j]->trackId()).pdgId()==11) {
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
	++nCellsSig;
      }
    }
    
    tntpArray[idx++] = nCells;
    tntpArray[idx++] = nCellsSig;
    tntpArray[idx++] = eDep;
    tntpArray[idx++] = eDepSig;
    
    _nCells=nCellsSig;

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
      double totEDepInWires(0.0), totEDepInWiresSig(0.0);
    
      for (size_t i=0; i<wSSteps; ++i) {
	// cout << "step " << i << ": looking at sense wire steps" << endl;
	StepPointMC const& whit = (*wShits)[i];
	SimParticleCollection::key_type trackId = whit.trackId();
	SimParticle const* sim = simParticles->getOrNull(trackId);
	if( !sim ) continue;
	wires.insert(whit.volumeId());
	totPathInWires+=whit.stepLength();
	totEDepInWires+=whit.eDep();
	*wirestxt << "Hit n. " << i+1 << " in the sense wire " << whit.volumeId()
		  << " made by " << sim->pdgId() << endl;
	if (sim->isPrimary() && sim->pdgId()==11 /*sim->fromGenerator()*/) {
	  wSStepsSig++;
	  sigwires.insert(whit.volumeId());
	  totPathInWiresSig+=whit.stepLength();
	  totEDepInWiresSig+=whit.eDep();
	}
      }
      tntpArray[idx++] = wires.size();
      tntpArray[idx++] = nSWHits;
      tntpArray[idx++] = totPathInWires;
      tntpArray[idx++] = totEDepInWires;
      tntpArray[idx++] = sigwires.size();
      tntpArray[idx++] = wSStepsSig;
      tntpArray[idx++] = totPathInWiresSig;
      tntpArray[idx++] = totEDepInWiresSig;
      //_tNtup->Fill(tntpArray);
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
    }
    
    
    if (nWlHits>0) {
      
      //      cout << "Some field wire has been hit!" << endl;
      
      size_t wlSteps = wllhits->size();
      
      //      cout << "Precisely they are " << wlSteps << endl;
      
      float wlStepsSig(0);
      
      set<unsigned long> walls, sigwalls, prevSigwalls;
      unsigned long old_wallId(0);
      bool noFirstWall(false);
      std::pair< set<unsigned long>::iterator, bool > insertInSW;
      double totPathInWalls(0.0), totPathInWallsSig(0.0);
      double totEDepInWalls(0.0), totEDepInWallsSig(0.0);
      size_t iHit=0;
      _nHit=0;

      for (size_t i=0; i<wlSteps; ++i) {
	//std::cout << "step " << i << ": looking at wall steps" << std::endl;
	StepPointMC const& whit = (*wllhits)[i];
	SimParticleCollection::key_type trackId = whit.trackId();
	SimParticle const* sim = simParticles->getOrNull(trackId);
	if( !sim ) continue;
	walls.insert(whit.volumeId());
	totPathInWalls+=whit.stepLength();
	totEDepInWalls+=whit.eDep();
	*wirestxt << "Hit n. " << i+1 << " in the wall " << whit.volumeId()
		  << " made by " << sim->pdgId() << endl;
	if (sim->isPrimary() && sim->pdgId()==11 /*sim->fromGenerator()*/) {
	  //std::cout<<"sgn el strp in wall: "<<whit.volumeId()<<" with eloss "<<whit.eDep()<<" dx "<<whit.stepLength()<<std::endl;
	  wlStepsSig++;
	  insertInSW = sigwalls.insert(whit.volumeId());
	  //std::cout<<"volID "<<whit.volumeId()<<" is already inserted? "<<insertInSW.second<<" is the first time? "<<(prevSigwalls.find(*insertInSW.first)==prevSigwalls.end())<<std::endl;
	  if (iHit<1000) {
	          if (insertInSW.second){
	                  //_nHit=sigwalls.size();
	                  //iHit=_nHit-1;
	                  iHit=_nHit;
	                  _wallEloss[iHit]=whit.eDep();
	                  _wallPath[iHit]=whit.stepLength();
	                  ++_nHit;
	                  if (noFirstWall) {
	                          prevSigwalls.insert(old_wallId);
	                  }
	                  old_wallId = *insertInSW.first;
	                  noFirstWall=true;
	          } else if ( prevSigwalls.find(*insertInSW.first)==prevSigwalls.end() ) {
	                  iHit = std::distance(sigwalls.begin(),insertInSW.first);
	                  _wallEloss[iHit]+=whit.eDep();
	                  _wallPath[iHit]+=whit.stepLength();
	          }
	  }
	  totPathInWallsSig+=whit.stepLength();
	  totEDepInWallsSig+=whit.eDep();
	}
      }
      tntpArray[idx++] = walls.size();
      tntpArray[idx++] = nWlHits;
      tntpArray[idx++] = totPathInWalls;
      tntpArray[idx++] = totEDepInWalls;
      tntpArray[idx++] = sigwalls.size();
      tntpArray[idx++] = wlStepsSig;
      tntpArray[idx++] = totPathInWallsSig;
      tntpArray[idx++] = totEDepInWallsSig;
      //_tNtup->Fill(tntpArray);
    } else {
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;
      tntpArray[idx++] = 0;

      _nHit=0;
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
	  if (id == VirtualDetectorId::TT_InSurf||
              id==VirtualDetectorId::TT_FrontPA||
              id==VirtualDetectorId::TT_Back) {
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
    _tSingleWallData->Fill();
    _tSingleCellData->Fill();
    
    delete [] tntpArray;
    
  } // end of doTracker
  
}

using mu2e::TSig;
DEFINE_ART_MODULE(TSig);

