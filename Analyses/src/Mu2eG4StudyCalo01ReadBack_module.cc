//
// Plugin to read/analyze g4study output
//
//  $Id: Mu2eG4StudyCalo01ReadBack_module.cc,v 1.5 2013/10/21 20:44:04 genser Exp $
//  $Author: genser $
//  $Date: 2013/10/21 20:44:04 $
//
// Original author KLG based on Mu2eG4StudyReadBack_module
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  typedef struct {

    Int_t run;
    Int_t evt;
    Int_t trk;
    Int_t endvol;

    Int_t pdg;
    Float_t time;
    Float_t gtime;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t p;
    Int_t code;

    Bool_t isstop;
    Float_t tstop;
    Float_t gtstop;
    Float_t xstop;
    Float_t ystop;
    Float_t zstop;
    Float_t pxstop;
    Float_t pystop;
    Float_t pzstop;
    Float_t pstop;
    Float_t pxprestop;
    Float_t pyprestop;
    Float_t pzprestop;
    Float_t pprestop;
    Int_t codestop;

    Int_t parent_id;
    Int_t parent_pdg;
    Float_t parent_x;
    Float_t parent_y;
    Float_t parent_z;
    Float_t parent_px;
    Float_t parent_py;
    Float_t parent_pz;
    Float_t parent_p;

  } NtPartData;

  typedef struct {

    Int_t run;
    Int_t evt;
    Int_t trk;
    Int_t vol;

    Int_t pdg;
    Float_t time;
    Float_t gtime;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t p;

    Float_t ke;
    Float_t tedep;
    Float_t niedep;

    Int_t endcode;



  } NtStepData;

  class Mu2eG4StudyCalo01ReadBack : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit Mu2eG4StudyCalo01ReadBack(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _stepperStepPoints(pset.get<string>("stepperStepPoints","stepper")),
      _tvdStepPoints(pset.get<string>("tvdStepPoints","timeVD")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _nttvd(0), _tpart(0), _tstep(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _timeCut(pset.get<double>("timeCut",0.0)),
      _stopped_only(pset.get<bool>("saveStopped",false)),
      _add_proper_time(pset.get<bool>("addProperTime",false))
    {

      nt = new float[128];

    }

    virtual ~Mu2eG4StudyCalo01ReadBack() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Names of the stepper and time virtual detectors StepPoint collections
    std::string _stepperStepPoints;
    std::string _tvdStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntstepper;
    TNtuple* _nttvd;
    TTree*   _tpart;
    TTree*   _tstep;

    float *nt; // Ntuple buffer
    NtPartData ttp;  // Buffer to fill particle tree
    NtStepData tstp; // Buffer to fill stepper tree

    // Pointers to the physical volumes we are interested in
    map<int,int> vid_stop;
    // we look at all of them for now

    // List of particles of interest for the particles ntuple
    //    set<int> pdg_save;
    // we save them all for now

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Save in the particles ntuple only those particles, which die 
    // after this time (in ns)
    double _timeCut;

    // Save only stopped particles in the particles ntuple
    bool _stopped_only;

    // Should we add together proper time for the whole decay chain
    bool _add_proper_time;

  };

  void Mu2eG4StudyCalo01ReadBack::beginJob(){

    vid_stop.clear();

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntstepper = tfs->make<TNtuple>( "ntstepper", "stepper ntuple",
                                     "evt:trk:sid:pdg:time:x:y:z:px:py:pz:gtime:code:run:ke:tedep:niedep");
    //                                0   1   2   3   4    5 6 7 8  9  10 11    12   13  14 15    16

    _nttvd = tfs->make<TNtuple>( "nttvd", "Time Virtual Detectors ntuple",
                                 "evt:trk:sid:pdg:time:x:y:z:px:py:pz:gtime:code:run:ke:tedep:niedep");
    //                            0   1   2   3   4    5 6 7 8  9  10 11    12   13  14 15    16

    _tpart = tfs->make<TTree>("tpart", "Particle tree");

    _tpart->Branch("run",        &ttp.run,        "run/I");
    _tpart->Branch("evt",        &ttp.evt,        "evt/I");
    _tpart->Branch("trk",        &ttp.trk,        "trk/I");
    _tpart->Branch("endvol",     &ttp.endvol,     "endvol/I");
    _tpart->Branch("pdg",        &ttp.pdg,        "pdg/I");
    _tpart->Branch("time",       &ttp.time,       "time/F");
    _tpart->Branch("gtime",      &ttp.gtime,      "gtime/F");
    _tpart->Branch("x",          &ttp.x,          "x/F");
    _tpart->Branch("y",          &ttp.y,          "y/F");
    _tpart->Branch("z",          &ttp.z,          "z/F");
    _tpart->Branch("px",         &ttp.px,         "px/F");
    _tpart->Branch("py",         &ttp.py,         "py/F");
    _tpart->Branch("pz",         &ttp.pz,         "pz/F");
    _tpart->Branch("p",          &ttp.p,          "p/F");
    _tpart->Branch("code",       &ttp.code,       "code/I");
    _tpart->Branch("isstop",     &ttp.isstop,     "isstop/O");
    _tpart->Branch("tstop",      &ttp.tstop,      "tstop/F");
    _tpart->Branch("gtstop",     &ttp.gtstop,     "gtstop/F");
    _tpart->Branch("xstop",      &ttp.xstop,      "xstop/F");
    _tpart->Branch("ystop",      &ttp.ystop,      "ystop/F");
    _tpart->Branch("zstop",      &ttp.zstop,      "zstop/F");
    _tpart->Branch("pxprestop",  &ttp.pxprestop,  "pxprestop/F");
    _tpart->Branch("pyprestop",  &ttp.pyprestop,  "pyprestop/F");
    _tpart->Branch("pzprestop",  &ttp.pzprestop,  "pzprestop/F");
    _tpart->Branch("pprestop",   &ttp.pprestop,   "pprestop/F");
    _tpart->Branch("pxstop",     &ttp.pxstop,     "pxstop/F");
    _tpart->Branch("pystop",     &ttp.pystop,     "pystop/F");
    _tpart->Branch("pzstop",     &ttp.pzstop,     "pzstop/F");
    _tpart->Branch("pstop",      &ttp.pstop,      "pstop/F");
    _tpart->Branch("codestop",   &ttp.codestop,   "codestop/I");
    _tpart->Branch("parent_id",  &ttp.parent_id,  "parent_id/I");
    _tpart->Branch("parent_pdg", &ttp.parent_pdg, "parent_pdg/I");
    _tpart->Branch("parent_x",   &ttp.parent_x,   "parent_x/F");
    _tpart->Branch("parent_y",   &ttp.parent_y,   "parent_y/F");
    _tpart->Branch("parent_z",   &ttp.parent_z,   "parent_z/F");
    _tpart->Branch("parent_px",  &ttp.parent_px,  "parent_px/F");
    _tpart->Branch("parent_py",  &ttp.parent_py,  "parent_py/F");
    _tpart->Branch("parent_pz",  &ttp.parent_pz,  "parent_pz/F");
    _tpart->Branch("parent_p",   &ttp.parent_p,   "parent_p/F");

    _tstep = tfs->make<TTree>("tstep", "Stepper tree");

    _tstep->Branch("run",        &tstp.run,        "run/I");
    _tstep->Branch("evt",        &tstp.evt,        "evt/I");
    _tstep->Branch("trk",        &tstp.trk,        "trk/I");
    _tstep->Branch("vol",        &tstp.vol,        "vol/I");
    _tstep->Branch("pdg",        &tstp.pdg,        "pdg/I");
    _tstep->Branch("time",       &tstp.time,       "time/F");
    _tstep->Branch("gtime",      &tstp.gtime,      "gtime/F");
    _tstep->Branch("x",          &tstp.x,          "x/F");
    _tstep->Branch("y",          &tstp.y,          "y/F");
    _tstep->Branch("z",          &tstp.z,          "z/F");
    _tstep->Branch("px",         &tstp.px,         "px/F");
    _tstep->Branch("py",         &tstp.py,         "py/F");
    _tstep->Branch("pz",         &tstp.pz,         "pz/F");
    _tstep->Branch("p",          &tstp.p,          "p/F");
    _tstep->Branch("ke",         &tstp.ke,         "ke/F");
    _tstep->Branch("tedep",      &tstp.tedep,      "tedep/F");
    _tstep->Branch("niedep",     &tstp.niedep,     "niedep/F");
    _tstep->Branch("endcode",    &tstp.endcode,    "endcode/I");

  }

  void Mu2eG4StudyCalo01ReadBack::beginRun(art::Run const& run){
  }

  void Mu2eG4StudyCalo01ReadBack::analyze(const art::Event& event) {

    ++_nAnalyzed;

    GlobalConstantsHandle<ParticleDataTable> pdt;
    ParticleDataTable const & pdt_ = *pdt;

    // print the pdt content

    static bool oneTime = true;

    if (oneTime) {

      oneTime = false;

      art::ServiceHandle<GeometryService> geom;
      SimpleConfig const& config  = geom->config();
      if (config.getBool("mu2e.printParticleDataTable",false)) {

        cout << __func__ 
             << " pdt size : "
             << pdt_.size() 
             << endl;
      
        for ( ParticleDataTable::const_iterator pdti=pdt_.begin(), e=pdt_.end(); 
              pdti!=e; ++pdti ) {
      
          cout << __func__ 
               << " pdt particle : "
               << pdti->first.pid()  
               << ", name: "          
               << pdt_.particle(pdti->first.pid()).ref().name()
               << ", PDTname: "          
               << pdt_.particle(pdti->first.pid()).ref().PDTname()
               << ", "
               << pdt_.particle(pdti->first.pid()).ref().mass()
               << ", "
               << pdt_.particle(pdti->first.pid()).ref().totalWidth()
               << ", "
               << pdt_.particle(pdti->first.pid()).ref().lifetime()
               << endl;

        }

      }

    } // end of oneTime

    // Ask the event to give us a "handle" to the requested points.
    art::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_stepperStepPoints,points);

    art::Handle<StepPointMCCollection> thits;
    event.getByLabel(_g4ModuleLabel,_tvdStepPoints,thits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());
    
    // Loop over all stepper points.
    if( points.isValid() ) for ( size_t i=0; i<points->size(); ++i ){

      // Alias, used for readability.
      StepPointMC const& point = (*points)[i];

      // Get the point information.

      CLHEP::Hep3Vector const& pos = point.position();
      CLHEP::Hep3Vector const& mom = point.momentum();

      // Get track info
      key_type trackId = point.trackId();
      int pdgId = 0;
      double mass(0.0);
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
      	  mass = pdt_.particle(pdgId).ref().mass();
	}
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId.asInt();
      nt[2]  = point.volumeId();
      nt[3]  = pdgId;
      nt[4]  = point.time();
      nt[5]  = pos.x();
      nt[6]  = pos.y();
      nt[7]  = pos.z();
      nt[8]  = mom.x();
      nt[9]  = mom.y();
      nt[10] = mom.z();
      nt[11] = point.properTime();
      nt[12] = point.endProcessCode();
      nt[13] = event.id().run();
      // compute kinetic energy: this is what Geant cuts on
      nt[14] = sqrt(mom.mag2()+mass*mass)-mass;
      nt[15] = point.totalEDep();
      nt[16] = point.nonIonizingEDep();

      _ntstepper->Fill(nt);
      if ( _nAnalyzed < _maxPrint){
        cout << "stepper point: "
             << event.id().run()   << " | "
             << event.id().event() << " | "
             << point.volumeId()   << " | "
             << point.trackId().asInt() << " | "
             << pdgId              << " , name: "  
             << pdt_.particle(pdgId).ref().name() << " , PDTname: "
             << pdt_.particle(pdgId).ref().PDTname() << " | "
             << point.time()       << " "
             << pos                << " "
             << mom.mag()          << " | "
             << point.totalEDep()  << " "
             << point.nonIonizingEDep()
             << endl;

      }

      // Fill the stepper tree (could "reuse" some of the nt/tree calc, or remove one or the other)

      tstp.run    = event.id().run();
      tstp.evt    = event.id().event();
      tstp.trk    = trackId.asInt();
      tstp.vol    = point.volumeId();
      tstp.pdg    = pdgId;
      tstp.time   = point.time();
      tstp.gtime  = point.properTime();
      tstp.x      = pos.x();
      tstp.y      = pos.y();
      tstp.z      = pos.z();  
      tstp.px     = mom.x(); 
      tstp.py     = mom.y(); 
      tstp.pz     = mom.z(); 
      tstp.p      = mom.mag();
      tstp.ke     = sqrt(mom.mag2()+mass*mass)-mass;
      tstp.tedep  = point.totalEDep();
      tstp.niedep = point.nonIonizingEDep();
      tstp.endcode= point.endProcessCode();

      _tstep->Fill();


    } // end loop over points.


    // Fill time virtual detectors ntuple

    // Loop over all time virtual detector hits.

    if( thits.isValid() ) for ( size_t i=0; i<thits->size(); ++i ){

      // Alias, used for readability.
      StepPointMC const& hit = (*thits)[i];

      // Get the hit information.

      int id = hit.volumeId();

      CLHEP::Hep3Vector const& pos = hit.position();
      CLHEP::Hep3Vector const& mom = hit.momentum();

      // Get track info
      key_type trackId = hit.trackId();
      int pdgId = 0;
      double mass(0.0);
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
	  mass = pdt_.particle(pdgId).ref().mass();
        }
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId.asInt();
      nt[2]  = id;
      nt[3]  = pdgId;
      nt[4]  = hit.time();
      nt[5]  = pos.x();
      nt[6]  = pos.y();
      nt[7]  = pos.z();
      nt[8]  = mom.x();
      nt[9]  = mom.y();
      nt[10] = mom.z();
      nt[11] = hit.properTime();
      nt[12] = hit.endProcessCode();
      nt[13] = event.id().run();
      nt[14] = sqrt(mom.mag2()+mass*mass)-mass;
      nt[15] = hit.totalEDep();
      nt[16] = hit.nonIonizingEDep();

      _nttvd->Fill(nt);

      if ( _nAnalyzed < _maxPrint){
        cout << "TVD hit:       "
             << event.id().run()   << " | "
             << event.id().event() << " | "
             << hit.volumeId()     << " | "
             << hit.trackId().asInt() << " | "
             << pdgId              << " , name: "  
             << pdt_.particle(pdgId).ref().name() << " , PDTname: "
             << pdt_.particle(pdgId).ref().PDTname() << " | "
             << hit.time()         << " "
             << pos                << " "
             << mom.mag()          << " | "
             << hit.totalEDep()    << " "
             << hit.nonIonizingEDep()
             << endl;
      }
    } // end loop over hits.

    // Fill tracks tree
    // with all sim partciles
    //    if( haveSimPart && pdg_save.size()>0 ) {
    if( haveSimPart ) {

      // Go through SimParticle container and analyze one particle at a time
      for ( SimParticleCollection::const_iterator isp=simParticles->begin();
            isp!=simParticles->end(); ++isp ){
        SimParticle const& sim = isp->second;

        // If particle PDG ID is not in the list - skip it
        // if( pdg_save.find(sim.pdgId()) == pdg_save.end() ) continue;
        // we save them all for now

        // Save SimParticle header info
        ttp.run = event.id().run();      // run_id
        ttp.evt = event.id().event();    // event_id
        ttp.trk = sim.id().asInt();      // track_id
        ttp.endvol = sim.endVolumeIndex();  // track last volume

        ttp.pdg = sim.pdgId();           // PDG id

	// Calculate parent proper time
	double gtime_parent = 0.0;
	if( _add_proper_time ) {
	  SimParticle const* sim_parent = &sim;
	  while( sim_parent && sim_parent->hasParent() ) {
	    sim_parent = simParticles->getOrNull(sim_parent->parentId());
	    if( sim_parent && sim_parent->pdgId()==ttp.pdg && sim.endDefined() ) {
	      gtime_parent += sim_parent->endProperTime();
	    }
	  }
	}

        // Save SimParticle other info
        ttp.time = sim.startGlobalTime(); // start time
        ttp.gtime = gtime_parent+sim.startProperTime(); // start time
        CLHEP::Hep3Vector const & pos_start = sim.startPosition();
        CLHEP::Hep3Vector const & mom_start = sim.startMomentum();
        ttp.x = pos_start.x();
        ttp.y = pos_start.y();
        ttp.z = pos_start.z();
        ttp.px = mom_start.x();
        ttp.py = mom_start.y();
        ttp.pz = mom_start.z();
        ttp.p  = mom_start.mag();
        ttp.code = sim.creationCode();

        // Check id of the volume where the particle died
        if( sim.endDefined() ) {

          ttp.isstop = ( vid_stop.find(sim.endVolumeIndex()) != vid_stop.end() );
          ttp.tstop = sim.endGlobalTime();
          ttp.gtstop = gtime_parent+sim.endProperTime();
          CLHEP::Hep3Vector const & pos_end = sim.endPosition();
          CLHEP::Hep3Vector const & mom_end = sim.endMomentum();

          ttp.xstop  = pos_end.x();
          ttp.ystop  = pos_end.y();
          ttp.zstop  = pos_end.z();

          // calculation of the prestep info is more involved...
          // assume the step points are sorted by time by construction, get the last one
          size_t thei(-1);
          size_t trackiid = sim.id().asInt();     
          for ( size_t i=points->size()-1; i!=0; --i ){
            if ( trackiid == ((*points)[i]).trackId().asInt() ) {
              thei = i;
              break;
            }
          }
          // we may need to protect this
          CLHEP::Hep3Vector const& mom_preend = (*points)[thei].momentum();

          ttp.pxprestop = mom_preend.x();
          ttp.pyprestop = mom_preend.y();
          ttp.pzprestop = mom_preend.z();
          ttp.pprestop  = mom_preend.mag();

          ttp.pxstop = mom_end.x();
          ttp.pystop = mom_end.y();
          ttp.pzstop = mom_end.z();
          ttp.pstop  = mom_end.mag();
          ttp.codestop = sim.stoppingCode();
        } else {
          ttp.isstop = false;
          ttp.tstop  = 0.0;
          ttp.gtstop = 0.0;
          ttp.xstop  = 0.0;
          ttp.ystop  = 0.0;
          ttp.zstop  = 0.0;

          ttp.pxprestop = 0.0;
          ttp.pyprestop = 0.0;
          ttp.pzprestop = 0.0;
          ttp.pprestop  = 0.0;

          ttp.pxstop    = 0.0;
          ttp.pystop    = 0.0;
          ttp.pzstop    = 0.0;
          ttp.pstop     = 0.0;

          ttp.codestop = 0;
        }

        // Parent info
        SimParticle const* sim_parent = 0;
        if( sim.hasParent() ) {
          ttp.parent_id = sim.parentId().asInt();
          sim_parent = simParticles->getOrNull(sim.parentId());
        } else {
          ttp.parent_id = -1;
        }
        if( sim_parent ) {
          ttp.parent_pdg = sim_parent->pdgId();
          CLHEP::Hep3Vector const & pos_parent = sim_parent->startPosition();
          CLHEP::Hep3Vector const & mom_parent = sim_parent->startMomentum();
          ttp.parent_x = pos_parent.x();
          ttp.parent_y = pos_parent.y();
          ttp.parent_z = pos_parent.z();
          ttp.parent_px = mom_parent.x();
          ttp.parent_py = mom_parent.y();
          ttp.parent_pz = mom_parent.z();
          ttp.parent_p  = mom_parent.mag();
        } else {
          ttp.parent_pdg = 0;
          ttp.parent_x = 0;
          ttp.parent_y = 0;
          ttp.parent_z = 0;
          ttp.parent_px = 0;
          ttp.parent_py = 0;
          ttp.parent_pz = 0;
          ttp.parent_p = 0;
        }

        // Keep only stopped particles
        if( _stopped_only && !ttp.isstop ) continue;

        // Keep only those particles, which die late enough
        if( _timeCut>0.0 && ttp.tstop<_timeCut ) continue;

        _tpart->Fill();

      }
    }

  }

}  // end namespace mu2e

using mu2e::Mu2eG4StudyCalo01ReadBack;
DEFINE_ART_MODULE(Mu2eG4StudyCalo01ReadBack);
