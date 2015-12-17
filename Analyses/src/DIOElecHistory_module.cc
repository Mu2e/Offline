//
// A module to follow the conversion electron in the events
//
// $Id: DIOElecHistory_module.cc,v 1.7 2013/10/21 20:44:04 genser Exp $
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
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "TFile.h"
#include "TTree.h"
#include "Mu2eUtilities/inc/BkgElecUtilities.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {


  class DIOElecHistory : public art::EDAnalyzer {
  public:
    explicit DIOElecHistory(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloROlabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrylabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _vdStepPoints(pset.get<std::string>("vdStepPoints", "virtualdetector")),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _tNtup(0),
      _nBadG4Status(0),
      _nOverflow(0),
      _nKilled(0),
      _totalcputime(0),
      _totalrealtime(0)
    {
    }
    virtual ~DIOElecHistory() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    void doEvent(art::Event const& evt);


    // Diagnostic level
    int _diagLevel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    std::string _caloROlabel, _caloCrylabel, _vdStepPoints;

    double _minimumEnergy; //minimum energy deposition of hits

    TTree* _tNtup;

    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;

    Float_t B_evt, B_run; // --2--
    Float_t B_gentime, B_genx, B_geny, B_genz, B_gene, B_genp, B_gencosth, B_genphi, B_genfoil; // --9--
    Float_t B_totedep, B_hitx[10000], B_hity[10000], B_hitz[10000], B_hitp[10000], B_hitcosth[10000], B_hitphi[10000], B_hittime[10000]; // --8--
    Int_t B_nhit; // --1--
    // Float_t B_ncrys, B_totcryedep; 
    Float_t B_cry1x, B_cry1y, B_cry1z, B_cry1p, B_cry1costh, B_cry1phi, B_cry1time; // --9--
    Float_t B_deadx, B_deady, B_deadz, B_deadtime, B_deadvol; // --5--
    Float_t B_vdx, B_vdy, B_vdz, B_vdp, B_vdcosth, B_vdphi, B_vdtime; // --7--
    Int_t   B_ndau, B_isfirst[1000]; // --2--
    Float_t B_daux[1000], B_dauy[1000], B_dauz[1000], B_daup[1000], B_daucosth[1000], B_dauphi[1000], B_dautime[1000], B_daupdgid[1000]; //--8--               
    Int_t B_nstraw;
    Float_t B_strawidx[10000], B_strawlay[10000], B_strawplane[10000], B_strawsec[10000], B_strawx[10000], B_strawy[10000], B_strawz[10000];
  };


  void DIOElecHistory::beginJob( ) {
  }

  void DIOElecHistory::analyze(art::Event const& evt ) {

    static int ncalls(0);
    ++ncalls;

    //   cout << "ncalls = " << ncalls << endl;
    art::Handle<StatusG4> g4StatusHandle;
    evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting DIOElecHistory::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    if (g4Status.overflowSimParticles()) {
      ++_nOverflow;
      mf::LogError("G4")
        << "Aborting DIOElecHistory::analyze due to overflow of particles\n"
        << g4Status;
      return;
    }

    if (g4Status.nKilledStepLimit() > 0) {
      ++_nKilled;
      mf::LogError("G4")
        << "Aborting DIOElecHistory::analyze due to nkilledStepLimit reached\n"
        << g4Status;
      //      return;
    }

    _totalcputime += g4Status.cpuTime();
    _totalrealtime += g4Status.realTime();

    if (ncalls == 1) {

      art::ServiceHandle<art::TFileService> tfs;

      _tNtup        = tfs->make<TTree>("ConvElec", "ConvElec Info");

      _tNtup->Branch("evt", &B_evt, "evt/F");
      _tNtup->Branch("run", &B_run, "run/F");
      _tNtup->Branch("gentime", &B_gentime, "gentime/F");
      _tNtup->Branch("genx", &B_genx, "genx/F");
      _tNtup->Branch("geny", &B_geny, "geny/F");
      _tNtup->Branch("genz", &B_genz, "genz/F");
      _tNtup->Branch("gene", &B_gene, "gene/F");
      _tNtup->Branch("genp", &B_genp, "genp/F");
      _tNtup->Branch("gencosth", &B_gencosth, "gencosth/F");
      _tNtup->Branch("genphi", &B_genphi, "genphi/F");
      _tNtup->Branch("genfoil", &B_genfoil, "genfoil/F");
      _tNtup->Branch("nhit", &B_nhit, "nhit/I");
      _tNtup->Branch("totedeptr", &B_totedep, "totedeptr/F");
      _tNtup->Branch("hitx[nhit]", B_hitx, "hitx[nhit]/F");
      _tNtup->Branch("hity[nhit]", B_hity, "hity[nhit]/F");
      _tNtup->Branch("hitz[nhit]", B_hitz, "hitz[nhit]/F");
      _tNtup->Branch("hitp[nhit]", B_hitp, "hitp[nhit]/F");
      _tNtup->Branch("hitcosth[nhit]", B_hitcosth, "hitcosth[nhit]/F");
      _tNtup->Branch("hitphi[nhit]", B_hitphi, "hitphi[nhit]/F");
      _tNtup->Branch("hittime[nhit]", B_hittime, "hittime[nhit]/F");
      _tNtup->Branch("isfirst[nhit]", B_isfirst, "isfirst[nhit]/I");

      _tNtup->Branch("nstraw", &B_nstraw, "nstraw/I");
      _tNtup->Branch("strawx[nstraw]", B_strawx, "strawx[nstraw]/F");
      _tNtup->Branch("strawy[nstraw]", B_strawy, "strawy[nstraw]/F");
      _tNtup->Branch("strawz[nstraw]", B_strawz, "strawz[nstraw]/F");
      _tNtup->Branch("strawidx[nstraw]", B_strawidx, "strawidx[nstraw]/F");
      _tNtup->Branch("strawlay[nstraw]", B_strawlay, "strawlay[nstraw]/F");
      _tNtup->Branch("strawplane[nstraw]", B_strawplane, "strawplane[nstraw]/F");
      _tNtup->Branch("strawsec[nstraw]", B_strawsec, "strawsec[nstraw]/F");

      //      _tNtup->Branch("ncrys", &B_ncrys, "ncrys/F");
      //   _tNtup->Branch("totcryedep", &B_totcryedep, "totcryedep/F");
      _tNtup->Branch("cry1x", &B_cry1x, "cry1x/F");
      _tNtup->Branch("cry1y", &B_cry1y, "cry1y/F");
      _tNtup->Branch("cry1z", &B_cry1z, "cry1z/F");
      _tNtup->Branch("cry1p", &B_cry1p, "cry1p/F");
      _tNtup->Branch("cry1costh", &B_cry1costh, "cry1costh/F");
      _tNtup->Branch("cry1phi", &B_cry1phi, "cry1phi/F");
      _tNtup->Branch("cry1time", &B_cry1time, "cry1time/F");
      _tNtup->Branch("deadx", &B_deadx, "deadx/F");
      _tNtup->Branch("deady", &B_deady, "deady/F");
      _tNtup->Branch("deadz", &B_deadz, "deadz/F");
      _tNtup->Branch("deadtime", &B_deadtime, "deadtime/F");
      _tNtup->Branch("deadvol", &B_deadvol, "deadvol/F");
      _tNtup->Branch("vdx", &B_vdx, "vdx/F");
      _tNtup->Branch("vdy", &B_vdy, "vdy/F");
      _tNtup->Branch("vdz", &B_vdz, "vdz/F");
      _tNtup->Branch("vdp", &B_vdp, "vdp/F");
      _tNtup->Branch("vdcosth", &B_vdcosth, "vdcosth/F");
      _tNtup->Branch("vdphi", &B_vdphi, "vdphi/F");
      _tNtup->Branch("vdtime", &B_vdtime, "vdtime/F");
      _tNtup->Branch("ndau", &B_ndau, "ndau/I");
      _tNtup->Branch("daux[ndau]", B_daux, "daux[ndau]/F");
      _tNtup->Branch("dauy[ndau]", B_dauy, "dauy[ndau]/F");
      _tNtup->Branch("dauz[ndau]", B_dauz, "dauz[ndau]/F");
      _tNtup->Branch("daup[ndau]", B_daup, "daup[ndau]/F");
      _tNtup->Branch("daucosth[ndau]", B_daucosth, "daucosth[ndau]/F");
      _tNtup->Branch("dauphi[ndau]", B_dauphi, "dauphi[ndau]/F");
      _tNtup->Branch("dautime[ndau]", B_dautime, "dautime[ndau]/F");
      _tNtup->Branch("daupdgid[ndau]", B_daupdgid, "daupdgid[ndau]/F");

    }
    
    doEvent(evt);
    
    
  } // end of analyze

    void DIOElecHistory::endJob() {
    cout << "DIOElecHistory::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
	 << "\nDIOElecHistory::endJob Number of overflow events "
         << "due to too many particles in G4: "
         << _nOverflow
	 << "\nDIOElecHistory::endJob Number of events with killed particles "
         << "due to too many steps in G4: "
         << _nKilled
	 << "\nDIOElecHistory::endJob total CpuTime "
         << _totalcputime
	 << "\nDIOElecHistory::endJob total RealTime "
         << _totalrealtime
         << endl;
  }


  void DIOElecHistory::doEvent(art::Event const& evt) {

    BkgElecUtilities Bkgt(evt,_generatorModuleLabel,
                           _g4ModuleLabel, _trackerStepPoints,
			   _caloROlabel);

 
    if (!Bkgt.hasStepPointMC()) return;

    B_evt = evt.id().event();
    B_run = evt.run();
    
    const GenParticle& genBkg = Bkgt.genBkgElec();

    B_gentime    = genBkg.time();
    B_genx       = genBkg.position().x();
    B_geny       = genBkg.position().y();
    B_genz       = genBkg.position().z();
    B_gene       = genBkg.momentum().e();
    B_genp       = genBkg.momentum().vect().mag();
    B_gencosth   = genBkg.momentum().cosTheta();
    B_genphi     = genBkg.momentum().phi();

    GeomHandle<StoppingTarget> target;
    float nfoil = 0;
    for (int i=0; i<target->nFoils(); ++i) {
      TargetFoil const& foil = target->foil(i);
      if (genBkg.position().z() >= foil.centerInMu2e().z()-foil.halfThickness() &&
          genBkg.position().z() <= foil.centerInMu2e().z()+foil.halfThickness() )
        nfoil = i;
    }
    B_genfoil  = nfoil;
    B_nhit     = Bkgt.hasStepPointMC();
    B_totedep  = Bkgt.totEDep();


    const vector<size_t> & hitsidx = Bkgt.bkgElecHitsIdx();

    art::Handle<StepPointMCCollection> hits;
    evt.getByLabel(_g4ModuleLabel, _trackerStepPoints, hits);

    size_t firstidx = 10000;
    double temptime = 100000;

    for (size_t i=0; i<hitsidx.size(); ++i) {

      StepPointMC const& hit = (*hits)[hitsidx[i]];
      
      B_hitx[i]     = hit.position().x();
      B_hity[i]     = hit.position().y();
      B_hitz[i]     = hit.position().z();
      B_hitp[i]     = hit.momentum().mag();
      B_hitcosth[i] = hit.momentum().cosTheta();
      B_hitphi[i]   = hit.momentum().phi();
      B_hittime[i]  = hit.time();
     
      if (hit.time() < temptime) {
	temptime = hit.time();
	firstidx = i;
      }      
    }

    for (size_t i=0; i<hitsidx.size(); ++i) {
      if (i==firstidx) {
	B_isfirst[i] = 1; 
      } else {
	B_isfirst[i] = 0;
      }
    }

    const Tracker& tracker = getTrackerOrThrow();

    const vector<StrawIndex> & strawvec = Bkgt.bkgElecStrawIdx();
    B_nstraw = strawvec.size();
    if (strawvec.size() > 0) {
      //      art::Handle<StrawHitCollection> pdataHandle;
      //  evt.getByLabel(_makerModuleLabel,pdataHandle);
      //  StrawHitCollection const* sthits = pdataHandle.product();
      for (size_t i=0; i<strawvec.size(); ++i) {
	Straw stw = tracker.getStraw(strawvec[i]);
	B_strawidx[i] = stw.id().getStraw();
	B_strawlay[i] = stw.id().getLayerId().getLayer();
	B_strawplane[i] = stw.id().getPlaneId();
	B_strawsec[i] = stw.id().getPanelId().getPanel();
	const CLHEP::Hep3Vector stMidPoint3 = stw.getMidPoint();
	B_strawx[i] = stMidPoint3.getX();
	B_strawy[i] = stMidPoint3.getY();
	B_strawz[i] = stMidPoint3.getZ();
      }
    } 

    if (Bkgt.gotCaloHit()) {
      StepPointMC const & calo1 = Bkgt.firstCaloHit();
      
      B_cry1x     = calo1.position().x() + 3904.;
      B_cry1y     = calo1.position().y();
      B_cry1z     = calo1.position().z() - 10200.;
      B_cry1p     = calo1.momentum().mag();
      B_cry1costh = calo1.momentum().cosTheta();
      B_cry1phi   = calo1.momentum().phi();
      B_cry1time  = calo1.time();
    } else {
      cout << "Calorimeter hit not found in the module!! " << endl;
      B_cry1x     = 0;
      B_cry1y     = 0;
      B_cry1z     = 0;
      B_cry1p     = 0;
      B_cry1costh = 0;
      B_cry1phi   = 0;
      B_cry1time  = 0;
    }

    const SimParticle& simBKG = Bkgt.simBkgElec(); 

    B_deadx  = simBKG.endPosition().x() + 3904.;
    B_deady  = simBKG.endPosition().y();
    B_deadz  = simBKG.endPosition().z() - 10200.;
    B_deadtime  = simBKG.endGlobalTime();
    B_deadvol  = simBKG.endVolumeIndex();


    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;
    art::Handle<StepPointMCCollection> vdhits;
    evt.getByLabel(_g4ModuleLabel,_vdStepPoints,vdhits);
    if (!vdhits.isValid()) return; 
    double time = 100000;
    size_t vdindex = 0;
    
    for (size_t i=0; i<vdhits->size(); ++i) {
      
      const StepPointMC& hit = (*vdhits)[i];
      
      int vdid = hit.volumeId();
      if (vdid != 11 && vdid != 12) continue;
      if (hit.trackId() == simBKG.id()) {
	if (hit.time() < time) time = hit.time();
	vdindex = i;
      }
    }
    
    const StepPointMC& hit = (*vdhits)[vdindex];
    
    CLHEP::Hep3Vector vdpos = hit.position();
    CLHEP::Hep3Vector vdmom = hit.momentum();
    Float_t vdtime = hit.time();
    
    B_vdx = vdpos.x() + 3904.;
    B_vdy = vdpos.y();
    B_vdz = vdpos.z() - 10200.;
    B_vdp = vdmom.mag();
    B_vdcosth = vdmom.cosTheta();
    B_vdphi = vdmom.phi();
    B_vdtime = vdtime;

    vector<art::Ptr<SimParticle> > const & BKGdau = simBKG.daughters();

    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);
    
    B_ndau  = BKGdau.size();

    for (size_t dauidx = 0; dauidx < BKGdau.size(); ++dauidx) {
      
      SimParticle const& d = *BKGdau[dauidx];

      PhysicalVolumeInfo const& volInfo = volumes->at(d.startVolumeIndex());
    
      cout << "Daughter n. " << dauidx << ": " << endl;
      cout << "pdg " << d.pdgId() << "'\t pos: "
	   << d.startPosition().z() - 10200.
	   << "\t start volume: " << volInfo.name() << endl;

  
      B_daux[dauidx] = d.startPosition().x() + 3904.;
      B_dauy[dauidx] = d.startPosition().y();
      B_dauz[dauidx] = d.startPosition().z() - 10200.;
      B_daup[dauidx] = d.startMomentum().mag();
      B_daucosth[dauidx] = d.startMomentum().cosTheta();
      B_dauphi[dauidx] = d.startMomentum().phi();
      B_dautime[dauidx] = d.startGlobalTime();
      B_daupdgid[dauidx] = d.pdgId();
    }
    
  

    PhysicalVolumeInfo const& volInfo = volumes->at(simBKG.endVolumeIndex());

    cout << "Event " << evt.id().event() << " : \nConversion Electron "
         << "dead in " << simBKG.endPosition() << " in the volume "
         << volInfo.name() << '\t' << simBKG.endVolumeIndex() << " because of "
         << simBKG.stoppingCode() << endl;


    _tNtup->Fill();

  } // end of doEvent
  
}

using mu2e::DIOElecHistory;
DEFINE_ART_MODULE(DIOElecHistory);

