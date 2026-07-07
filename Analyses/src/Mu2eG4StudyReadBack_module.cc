//
// Plugin to read/analyze g4study output
//
//
// Original author KLG somewhat based on vd read back
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
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

  class Mu2eG4StudyReadBack : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit Mu2eG4StudyReadBack(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _stepperStepPoints(pset.get<string>("stepperStepPoints","stepper")),
      _tvdStepPoints(pset.get<string>("tvdStepPoints","timeVD")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _nttvd(0), _tpart(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _timeCut(pset.get<double>("timeCut",0.0)),
      _stopped_only(pset.get<bool>("saveStopped",false)),
      _add_proper_time(pset.get<bool>("addProperTime",false)),
      physVolInfoInput_(pset.get<std::string>("physVolInfoInput","g4run")),
      vols_(),
      verbosityLevel_(pset.get<int>("diagLevel", 0))
    {

      nt = new float[128];

    }

    virtual ~Mu2eG4StudyReadBack() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);
    virtual void beginSubRun(art::SubRun const&);

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

    float *nt; // Need this buffer to fill TTree ntvd
    NtPartData ttp; // Buffer to fill particle tree

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

    art::InputTag physVolInfoInput_;
    const PhysicalVolumeInfoMultiCollection *vols_;

    int verbosityLevel_;

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

  }

  void Mu2eG4StudyReadBack::beginRun(art::Run const& run){
  }

  void Mu2eG4StudyReadBack::beginSubRun(art::SubRun const& sr) {
    art::Handle<PhysicalVolumeInfoMultiCollection> volh;
    sr.getByLabel(physVolInfoInput_, volh);
    if( volh.isValid() ) {
      vols_ = &*volh;
        std::cout<<"PhysicalVolumeInfoMultiCollection dump begin"<<std::endl;
        for(const auto& i : *vols_) {
          std::cout<<"*********************************************************"<<std::endl;
          std::cout<<"Collection size = "<<i.size()<<std::endl;
          // register all volumes
          size_t ii=0;
          for(const auto& entry : i) {
            vid_stop[ii] = (entry.second).copyNo();
            std::cout<<entry.second<<" "<<ii<<std::endl;
            ++ii;
          }
        }
        std::cout<<"PhysicalVolumeInfoMultiCollection dump end"<<std::endl;
      }
  }

  void Mu2eG4StudyReadBack::analyze(const art::Event& event) {

    ++_nAnalyzed;

    GlobalConstantsHandle<ParticleDataList> pdt;
    ParticleDataList const & pdt_ = *pdt;

    // print the pdt content

    static bool oneTime = true;

    if (oneTime) {

      oneTime = false;

      art::ServiceHandle<GeometryService> geom;
      SimpleConfig const& config  = geom->config();
      if (config.getBool("mu2e.printParticleDataTable",false)) {
        pdt_.printTable();
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
          mass = pdt_.particle(pdgId).mass();
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
             << point.trackId().asInt() << " | "
             << pdgId              << " , name: "
             << pdt_.particle(pdgId).name() << " | "
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
          mass = pdt_.particle(pdgId).mass();
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
             << hit.trackId().asInt() << " | "
             << pdgId              << " , name: "
             << pdt_.particle(pdgId).name() << " | "
             << hit.time()         << " "
             << pos                << " "
             << mom.mag()
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
        CLHEP::Hep3Vector const mom_start = sim.startMomentum();
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
          CLHEP::Hep3Vector const mom_end = sim.endMomentum();

          ttp.xstop  = pos_end.x();
          ttp.ystop  = pos_end.y();
          ttp.zstop  = pos_end.z();

          ttp.pxstop = mom_end.x();
          ttp.pystop = mom_end.y();
          ttp.pzstop = mom_end.z();
          ttp.pstop  = mom_end.mag();
          ttp.codestop = sim.stoppingCode();

          // calculation of the prestep info is more involved...
          // assume the step points are sorted by time by construction, get the last one
          int thei(-1);
          size_t trackiid = sim.id().asInt();
          for ( int i=points->size()-1; i!=-1; --i ){
            // if ( _nAnalyzed < _maxPrint){
            //   cout << __func__ << " steppoint trackid = "
            //        << ((*points)[i]).trackId().asInt()
            //        << endl;
            // }
            if ( trackiid == ((*points)[i]).trackId().asInt() ) {
              thei = i;
              break;
            }
          }

          if ( thei >= 0 ) {

            CLHEP::Hep3Vector const& mom_preend = (*points)[thei].momentum();

            ttp.pxprestop = mom_preend.x();
            ttp.pyprestop = mom_preend.y();
            ttp.pzprestop = mom_preend.z();
            ttp.pprestop  = mom_preend.mag();

          } else {

            cout << __func__ << " WARNING thei = " << thei
                 << " did not find matching steppoint in event "
                 << event.id().event()
                 << " points->size() "
                 << points->size()
                 << " for trackiid "
                 << trackiid
                 << endl;

            ttp.pxprestop = 0.0;
            ttp.pyprestop = 0.0;
            ttp.pzprestop = 0.0;
            ttp.pprestop  = 0.0;

          }

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
          CLHEP::Hep3Vector const mom_parent = sim_parent->startMomentum();
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

using mu2e::Mu2eG4StudyReadBack;
DEFINE_ART_MODULE(Mu2eG4StudyReadBack)
