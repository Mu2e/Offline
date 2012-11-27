//
// Plugin to read/analyze g4study output
//
//  $Id: Mu2eG4StudyReadBack_module.cc,v 1.4 2012/11/27 23:00:59 genser Exp $
//  $Author: genser $
//  $Date: 2012/11/27 23:00:59 $
//
// Original author KLG somewhat based on vd read back
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  typedef struct {

    Int_t run;
    Int_t evt;
    Int_t trk;

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

  class Mu2eG4StudyReadBack : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit Mu2eG4StudyReadBack(fhicl::ParameterSet const& pset) :
      _stepperStepPoints(pset.get<string>("stepperStepPoints","stepper")),
      _tvdStepPoints(pset.get<string>("tvdStepPoints","timeVD")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _nttvd(0), _ntpart(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _timeCut(pset.get<double>("timeCut",0.0)),
      _stopped_only(pset.get<bool>("saveStopped",false)),
      _add_proper_time(pset.get<bool>("addProperTime",false))
    {

      nt = new float[128];

    }

    virtual ~Mu2eG4StudyReadBack() { }

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
    TTree* _ntpart;

    float *nt; // Need this buffer to fill TTree ntvd
    NtPartData ntp; // Buffer to fill particles ntuple

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

  void Mu2eG4StudyReadBack::beginJob(){

    vid_stop.clear();

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntstepper = tfs->make<TNtuple>( "ntstepper", "stepper ntuple",
                                     "evt:trk:sid:pdg:time:x:y:z:px:py:pz:gtime:code:run:ke");
    //                                0   1   2   3   4    5 6 7 8  9  10 11    12   13  14

    _nttvd = tfs->make<TNtuple>( "nttvd", "Time Virtual Detectors ntuple",
                                 "evt:trk:sid:pdg:time:x:y:z:px:py:pz:gtime:code:run:ke");
    //                            0   1   2   3   4    5 6 7 8  9  10 11    12   13  14

    _ntpart = tfs->make<TTree>("ntpart", "Particles ntuple");

    _ntpart->Branch("run",        &ntp.run,        "run/I");
    _ntpart->Branch("evt",        &ntp.evt,        "evt/I");
    _ntpart->Branch("trk",        &ntp.trk,        "trk/I");
    _ntpart->Branch("pdg",        &ntp.pdg,        "pdg/I");
    _ntpart->Branch("time",       &ntp.time,       "time/F");
    _ntpart->Branch("gtime",      &ntp.gtime,      "gtime/F");
    _ntpart->Branch("x",          &ntp.x,          "x/F");
    _ntpart->Branch("y",          &ntp.y,          "y/F");
    _ntpart->Branch("z",          &ntp.z,          "z/F");
    _ntpart->Branch("px",         &ntp.px,         "px/F");
    _ntpart->Branch("py",         &ntp.py,         "py/F");
    _ntpart->Branch("pz",         &ntp.pz,         "pz/F");
    _ntpart->Branch("p",          &ntp.p,          "p/F");
    _ntpart->Branch("code",       &ntp.code,       "code/I");
    _ntpart->Branch("isstop",     &ntp.isstop,     "isstop/O");
    _ntpart->Branch("tstop",      &ntp.tstop,      "tstop/F");
    _ntpart->Branch("gtstop",     &ntp.gtstop,     "gtstop/F");
    _ntpart->Branch("xstop",      &ntp.xstop,      "xstop/F");
    _ntpart->Branch("ystop",      &ntp.ystop,      "ystop/F");
    _ntpart->Branch("zstop",      &ntp.zstop,      "zstop/F");
    _ntpart->Branch("codestop",   &ntp.codestop,   "codestop/I");
    _ntpart->Branch("parent_id",  &ntp.parent_id,  "parent_id/I");
    _ntpart->Branch("parent_pdg", &ntp.parent_pdg, "parent_pdg/I");
    _ntpart->Branch("parent_x",   &ntp.parent_x,   "parent_x/F");
    _ntpart->Branch("parent_y",   &ntp.parent_y,   "parent_y/F");
    _ntpart->Branch("parent_z",   &ntp.parent_z,   "parent_z/F");
    _ntpart->Branch("parent_px",  &ntp.parent_px,  "parent_px/F");
    _ntpart->Branch("parent_py",  &ntp.parent_py,  "parent_py/F");
    _ntpart->Branch("parent_pz",  &ntp.parent_pz,  "parent_pz/F");
    _ntpart->Branch("parent_p",   &ntp.parent_p,   "parent_p/F");

  }

  void Mu2eG4StudyReadBack::beginRun(art::Run const& run){

    // Get pointers to the physical volumes we are interested; fixme not used for now
    art::Handle<PhysicalVolumeInfoCollection> physVolumes;
    run.getByLabel(_g4ModuleLabel, physVolumes);
    if( physVolumes.isValid() ) {

      // fixme get it from the config file

      // register all volumes
      for ( size_t i=0; i<physVolumes->size(); ++i ) {
        // if( (*physVolumes)[i].name() == "BoxInTheWorld" ) {
        vid_stop[i] = (*physVolumes)[i].copyNo();
        cout << "Mu2eG4StudyReadBack: register volume " << i << " = "
             << (*physVolumes)[i].name() << " " << (*physVolumes)[i].copyNo() << endl;
        //}
      }

    }

  }

  void Mu2eG4StudyReadBack::analyze(const art::Event& event) {

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
      const StepPointMC& point = (*points)[i];

      // Get the point information.

      const CLHEP::Hep3Vector& pos = point.position();
      const CLHEP::Hep3Vector& mom = point.momentum();

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

      _ntstepper->Fill(nt);
      if ( _nAnalyzed < _maxPrint){
        cout << "stepper point: "
             << event.id().run()   << " | "
             << event.id().event() << " | "
             << point.volumeId()   << " | "
             << pdgId              << " , name: "  
             << pdt_.particle(pdgId).ref().name() << " , PDTname: "
             << pdt_.particle(pdgId).ref().PDTname() << " | "
             << point.time()       << " "
             << pos                << " "
             << mom.mag()
             << endl;

      }

    } // end loop over points.


    // Fill time virtual detectors ntuple

    // Loop over all time virtual detector hits.

    if( thits.isValid() ) for ( size_t i=0; i<thits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*thits)[i];

      // Get the hit information.

      int id = hit.volumeId();

      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

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

      _nttvd->Fill(nt);

      if ( _nAnalyzed < _maxPrint){
        cout << "TVD hit:       "
             << event.id().run()   << " | "
             << event.id().event() << " | "
             << hit.volumeId()     << " | "
             << pdgId              << " , name: "  
             << pdt_.particle(pdgId).ref().name() << " , PDTname: "
             << pdt_.particle(pdgId).ref().PDTname() << " | "
             << hit.time()         << " "
             << pos                << " "
             << mom.mag()
             << endl;
      }
    } // end loop over hits.

    // Fill tracks ntuple
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
        ntp.run = event.id().run();      // run_id
        ntp.evt = event.id().event();    // event_id
        ntp.trk = sim.id().asInt();      // track_id
        ntp.pdg = sim.pdgId();           // PDG id

	// Calculate parent proper time
	double gtime_parent = 0.0;
	if( _add_proper_time ) {
	  SimParticle const* sim_parent = &sim;
	  while( sim_parent && sim_parent->hasParent() ) {
	    sim_parent = simParticles->getOrNull(sim_parent->parentId());
	    if( sim_parent && sim_parent->pdgId()==ntp.pdg && sim.endDefined() ) {
	      gtime_parent += sim_parent->endProperTime();
	    }
	  }
	}

        // Save SimParticle other info
        ntp.time = sim.startGlobalTime(); // start time
        ntp.gtime = gtime_parent+sim.startProperTime(); // start time
        CLHEP::Hep3Vector const & pos_start = sim.startPosition();
        CLHEP::Hep3Vector const & mom_start = sim.startMomentum();
        ntp.x = pos_start.x();
        ntp.y = pos_start.y();
        ntp.z = pos_start.z();
        ntp.px = mom_start.x();
        ntp.py = mom_start.y();
        ntp.pz = mom_start.z();
        ntp.p = mom_start.mag();
        ntp.code = sim.creationCode();

        // Check id of the volume where the particle died
        if( sim.endDefined() ) {
          if( vid_stop.find(sim.endVolumeIndex()) != vid_stop.end() ) {
            ntp.isstop = true;
          } else {
            ntp.isstop = false;
          }
          ntp.tstop = sim.endGlobalTime();
          ntp.gtstop = gtime_parent+sim.endProperTime();
          CLHEP::Hep3Vector const & pos_end = sim.endPosition();
          ntp.xstop = pos_end.x();
          ntp.ystop = pos_end.y();
          ntp.zstop = pos_end.z();
          ntp.codestop = sim.stoppingCode();
        } else {
          ntp.isstop = false;
          ntp.tstop = 0;
          ntp.gtstop = 0;
          ntp.xstop = 0;
          ntp.ystop = 0;
          ntp.zstop = 0;
          ntp.codestop = 0;
        }

        // Parent info
        SimParticle const* sim_parent = 0;
        if( sim.hasParent() ) {
          ntp.parent_id = sim.parentId().asInt();
          sim_parent = simParticles->getOrNull(sim.parentId());
        } else {
          ntp.parent_id = -1;
        }
        if( sim_parent ) {
          ntp.parent_pdg = sim_parent->pdgId();
          CLHEP::Hep3Vector const & pos_parent = sim_parent->startPosition();
          CLHEP::Hep3Vector const & mom_parent = sim_parent->startMomentum();
          ntp.parent_x = pos_parent.x();
          ntp.parent_y = pos_parent.y();
          ntp.parent_z = pos_parent.z();
          ntp.parent_px = mom_parent.x();
          ntp.parent_py = mom_parent.y();
          ntp.parent_pz = mom_parent.z();
          ntp.parent_p = mom_parent.mag();
        } else {
          ntp.parent_pdg = 0;
          ntp.parent_x = 0;
          ntp.parent_y = 0;
          ntp.parent_z = 0;
          ntp.parent_px = 0;
          ntp.parent_py = 0;
          ntp.parent_pz = 0;
          ntp.parent_p = 0;
        }

        // Keep only stopped particles
        if( _stopped_only && !ntp.isstop ) continue;

        // Keep only those particles, which die late enough
        if( _timeCut>0.1 && ntp.tstop<_timeCut ) continue;

        _ntpart->Fill();

      }
    }

  }

}  // end namespace mu2e

using mu2e::Mu2eG4StudyReadBack;
DEFINE_ART_MODULE(Mu2eG4StudyReadBack);
