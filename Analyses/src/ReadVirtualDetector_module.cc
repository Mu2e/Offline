//
// Plugin to read virtual detectors data and create ntuples
//
//  $Id: ReadVirtualDetector_module.cc,v 1.7 2011/07/17 01:21:58 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/07/17 01:21:58 $
//
// Original author Ivan Logashenko
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Persistency/Common/Handle.h"
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

  const unsigned int nvdet = 20;

  typedef struct {

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

    Int_t g4bl_evt;
    Int_t g4bl_trk;
    Float_t g4bl_weight;
    Float_t g4bl_time;

    Int_t parent_id;
    Int_t parent_pdg;
    Float_t parent_x;
    Float_t parent_y;
    Float_t parent_z;
    Float_t parent_px;
    Float_t parent_py;
    Float_t parent_pz;
    Float_t parent_p;

    Int_t nvd;
    Bool_t isvd[nvdet];
    Float_t tvd[nvdet];
    Float_t gtvd[nvdet];
    Float_t xvd[nvdet];
    Float_t yvd[nvdet];
    Float_t zvd[nvdet];
    Float_t pxvd[nvdet];
    Float_t pyvd[nvdet];
    Float_t pzvd[nvdet];
    Float_t pvd[nvdet];
    Float_t xlvd[nvdet];
    Float_t ylvd[nvdet];
    Float_t zlvd[nvdet];

  } NtPartData;

  class ReadVirtualDetector : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;
    typedef SimParticleCollection::key_type key_type;

    explicit ReadVirtualDetector(fhicl::ParameterSet const& pset) :
      _vdStepPoints(pset.get<string>("vdStepPoints","virtualdetector")),
      _tvdStepPoints(pset.get<string>("tvdStepPoints","timeVD")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _ntvd(0), _nttvd(0), _ntpart(0), _ntpart1(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _vd_required(pset.get<int>("requireVD",0)),
      _timeCut(pset.get<double>("timeCut",0.0)),
      _stopped_only(pset.get<bool>("saveStopped",false))
    {

      Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
      if( pdg_ids.size()>0 ) {
        cout << "ReadVirtualDetector: save following particle types in the ntuple: ";
        for( size_t i=0; i<pdg_ids.size(); ++i ) {
          pdg_save.insert(pdg_ids[i]);
          cout << pdg_ids[i] << ", ";
        }
        cout << endl;
      }

      Vint const & vd_ids = pset.get<Vint>("saveVD", Vint());
      if( vd_ids.size()>0 ) {
        cout << "ReadVirtualDetector: save data from the following virtual detectors: ";
        for( size_t i=0; i<vd_ids.size(); ++i ) {
          vd_save.insert(vd_ids[i]);
          cout << vd_ids[i] << ", ";
        }
        cout << endl;
      }

      nt = new float[1000];

    }

    virtual ~ReadVirtualDetector() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Name of the VD and TVD StepPoint collections
    std::string  _vdStepPoints;
    std::string _tvdStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntvd;
    TNtuple* _nttvd;
    TTree* _ntpart;
    TTree* _ntpart1;

    float *nt; // Need this buffer to fill TTree ntvd
    NtPartData ntp; // Buffer to fill particles ntuple

    // Pointers to the physical volumes we are interested in
    // -- stopping target
    map<int,int> vid_stop;

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // List of virtual detectors to be saved
    set<int> vd_save;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Virtual detector, which has to be crossed by particle before
    // it is saved in particles ntuple
    int _vd_required;

    // Save in the particles ntuple only those particles, which die 
    // after this time (in ns)
    double _timeCut;

    // Save only stopped particles in the particles ntuple
    bool _stopped_only;

  };

  void ReadVirtualDetector::beginJob(){

    vid_stop.clear();

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntvd = tfs->make<TNtuple>( "ntvd", "Virtual Detectors ntuple",
                                "evt:trk:sid:pdg:time:x:y:z:px:py:pz:"
                                "xl:yl:zl:pxl:pyl:pzl:gtime:"
                                "g4bl_weight:g4bl_time");

    _nttvd = tfs->make<TNtuple>( "nttvd", "Time Virtual Detectors ntuple",
                                 "evt:trk:sid:pdg:time:x:y:z:px:py:pz:"
                                 "gtime:code:g4bl_weight:g4bl_time");

    // Have to use TTree here, because one cannot use more than 100 variables in TNtuple

    _ntpart = tfs->make<TTree>("ntpart", "Particles ntuple");
    /*
      _ntpart->Branch("all",nt,
      "evt:trk:pdg:"
      "time:gtime:x:y:z:px:py:pz:"
      "isstop:tstop:gtstop:xstop:ystop:zstop:"
      "g4bl_evt:g4bl_trk:g4bl_weight:g4bl_time:"
      "parent_id:parent_pdg:"
      "parent_x:parent_y:parent_z:"
      "parent_px:parent_py:parent_pz:"
      "nvd:isvd[20]:"
      "tvd[20]:gtvd[20]:xvd[20]:yvd[20]:zvd[20]:"
      "pxvd[20]:pyvd[20]:pzvd[20]:"
      "xlvd[20]:ylvd[20]:zlvd[20]"
      );
    */

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
    _ntpart->Branch("g4bl_evt",   &ntp.g4bl_evt,   "g4bl_evt/I");
    _ntpart->Branch("g4bl_trk",   &ntp.g4bl_trk,   "g4bl_trk/I");
    _ntpart->Branch("g4bl_weight",&ntp.g4bl_weight,"g4bl_weight/F");
    _ntpart->Branch("g4bl_time",  &ntp.g4bl_time,  "g4bl_time/F");
    _ntpart->Branch("parent_id",  &ntp.parent_id,  "parent_id/I");
    _ntpart->Branch("parent_pdg", &ntp.parent_pdg, "parent_pdg/I");
    _ntpart->Branch("parent_x",   &ntp.parent_x,   "parent_x/F");
    _ntpart->Branch("parent_y",   &ntp.parent_y,   "parent_y/F");
    _ntpart->Branch("parent_z",   &ntp.parent_z,   "parent_z/F");
    _ntpart->Branch("parent_px",  &ntp.parent_px,  "parent_px/F");
    _ntpart->Branch("parent_py",  &ntp.parent_py,  "parent_py/F");
    _ntpart->Branch("parent_pz",  &ntp.parent_pz,  "parent_pz/F");
    _ntpart->Branch("parent_p",   &ntp.parent_p,   "parent_p/F");
    _ntpart->Branch("nvd",        &ntp.nvd,        "nvd/I");
    _ntpart->Branch("isvd",        ntp.isvd,       "isvd[nvd]/F");
    _ntpart->Branch("tvd",         ntp.tvd,        "tvd[nvd]/F");
    _ntpart->Branch("gtvd",        ntp.gtvd,       "gtvd[nvd]/F");
    _ntpart->Branch("xvd",         ntp.xvd,        "xvd[nvd]/F");
    _ntpart->Branch("yvd",         ntp.yvd,        "yvd[nvd]/F");
    _ntpart->Branch("zvd",         ntp.zvd,        "zvd[nvd]/F");
    _ntpart->Branch("pxvd",        ntp.pxvd,       "pxvd[nvd]/F");
    _ntpart->Branch("pyvd",        ntp.pyvd,       "pyvd[nvd]/F");
    _ntpart->Branch("pzvd",        ntp.pzvd,       "pzvd[nvd]/F");
    _ntpart->Branch("pvd",         ntp.pvd,        "pvd[nvd]/F");
    _ntpart->Branch("xlvd",        ntp.xlvd,       "xlvd[nvd]/F");
    _ntpart->Branch("ylvd",        ntp.ylvd,       "ylvd[nvd]/F");
    _ntpart->Branch("zlvd",        ntp.zlvd,       "zlvd[nvd]/F");

  }

  void ReadVirtualDetector::beginRun(art::Run const& run){

    // Get pointers to the physical volumes we are interested
    art::Handle<PhysicalVolumeInfoCollection> physVolumes;
    run.getByLabel(_g4ModuleLabel, physVolumes);
    if( physVolumes.isValid() ) {

      for ( size_t i=0; i<physVolumes->size(); ++i ) {
        if( (*physVolumes)[i].name() == "TargetFoil_" ) {
          vid_stop[i] = (*physVolumes)[i].copyNo();
          cout << "ReadVirtualDetector: register stopping target volume " << i << " = "
               << (*physVolumes)[i].name() << " " << (*physVolumes)[i].copyNo() << endl;
        }
      }

    }

  }

  void ReadVirtualDetector::analyze(const art::Event& event) {

    ++_nAnalyzed;

    // Access virtual detectors geometry information
    // If not virtual detectors are defined, skip the rest

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_vdStepPoints,hits);

    art::Handle<StepPointMCCollection> thits;
    event.getByLabel(_g4ModuleLabel,_tvdStepPoints,thits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    art::Handle<G4BeamlineInfoCollection> g4beamlineData;
    event.getByLabel(_generatorModuleLabel, g4beamlineData);
    bool haveG4BL = g4beamlineData.isValid();
    if ( haveG4BL ) haveG4BL = (g4beamlineData->size()==1);

    // Fill virtual detectors ntuple

    if( haveG4BL ) {
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      nt[18] = extra.weight();
      nt[19] = extra.time();
    } else {
      nt[18] = 0;
      nt[19] = 0;
    }

    // Loop over all VD hits.
    if( hits.isValid() ) for ( size_t i=0; i<hits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.

      int id = hit.volumeId();

      // If virtual detector id is not in the list - skip it
      if( vd_save.size()>0 && vd_save.find(id) == vd_save.end() ) continue;

      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
      CLHEP::Hep3Vector lmom = mom;
      if( vdg->getRotation(id)!=0 ) {
        lpos *= *(vdg->getRotation(id));
        lmom *= *(vdg->getRotation(id));
      }

      // Get track info
      key_type trackId = hit.trackId();
      int pdgId = 0;
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
        }
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId.asInt();
      nt[2]  = hit.volumeId();
      nt[3]  = pdgId;
      nt[4]  = hit.time();
      nt[5]  = pos.x();
      nt[6]  = pos.y();
      nt[7]  = pos.z();
      nt[8]  = mom.x();
      nt[9]  = mom.y();
      nt[10] = mom.z();
      nt[11] = lpos.x();
      nt[12] = lpos.y();
      nt[13] = lpos.z();
      nt[14] = lmom.x();
      nt[15] = lmom.y();
      nt[16] = lmom.z();
      nt[17] = hit.properTime();

      _ntvd->Fill(nt);

      if ( _nAnalyzed < _maxPrint){
        cout << "VD hit: "
             << event.id().event() << " | "
             << hit.volumeId()     << " "
             << pdgId              << " | "
             << hit.time()         << " "
             << lpos               << " "
             << mom.mag()
             << endl;

      }

    } // end loop over hits.


    // Fill time virtual detectors ntuple

    if( haveG4BL ) {
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      nt[13] = extra.weight();
      nt[14] = extra.time();
    } else {
      nt[13] = 0;
      nt[14] = 0;
    }

    // Loop over all time VD hits.

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
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
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

      _nttvd->Fill(nt);

      if ( _nAnalyzed < _maxPrint){
        cout << "TVD hit: "
             << event.id().event() << " | "
             << hit.volumeId()     << " "
             << pdgId              << " | "
             << hit.time()         << " "
             << pos                << " "
             << mom.mag()
             << endl;

      }

    } // end loop over hits.

    // Fill tracks ntuple
    if( haveSimPart && pdg_save.size()>0 ) {


      // Go through SimParticle container and analyze one particle at a time
      for ( SimParticleCollection::const_iterator isp=simParticles->begin();
            isp!=simParticles->end(); ++isp ){
        SimParticle const& sim = isp->second;

        // It particle PDG ID is not in the list - skip it
        if( pdg_save.find(sim.pdgId()) == pdg_save.end() ) continue;

        // Save SimParticle info
        ntp.evt = event.id().event();    // event_id
        ntp.trk = sim.id().asInt();      // track_id
        ntp.pdg = sim.pdgId();           // PDG id
        ntp.time = sim.startGlobalTime(); // start time
        ntp.gtime = sim.startProperTime(); // start time
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

        // Check id of the volume there particle dies
        if( sim.endDefined() ) {
          if( vid_stop.find(sim.endVolumeIndex()) != vid_stop.end() ) {
            ntp.isstop = true;
          } else {
            ntp.isstop = false;
          }
          ntp.tstop = sim.endGlobalTime();
          ntp.gtstop = sim.endProperTime();
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

        if( haveG4BL ) {
          G4BeamlineInfo const& extra = g4beamlineData->at(0);
          ntp.g4bl_evt    = extra.eventId();
          ntp.g4bl_trk    = extra.trackId();
          ntp.g4bl_weight = extra.weight();
          ntp.g4bl_time   = extra.time();
        } else {
          ntp.g4bl_evt    = 0;
          ntp.g4bl_trk    = 0;
          ntp.g4bl_weight = 0;
          ntp.g4bl_time   = 0;
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

        // Clear up VD data
        ntp.nvd = nvdet;
        for ( size_t i=0; i<nvdet; ++i ) {
          ntp.isvd[i]=false;
          ntp.tvd[i]=0;
          ntp.gtvd[i]=0;
          ntp.xvd[i]=0;
          ntp.yvd[i]=0;
          ntp.zvd[i]=0;
          ntp.pxvd[i]=0;
          ntp.pyvd[i]=0;
          ntp.pzvd[i]=0;
          ntp.pvd[i]=0;
          ntp.xlvd[i]=0;
          ntp.ylvd[i]=0;
          ntp.zlvd[i]=0;
        }

        // Loop over all virtual detectors and fill corresponding data
        for ( size_t i=0; i<hits->size(); ++i ){

          // Alias, used for readability.
          const StepPointMC& hit = (*hits)[i];

          // Only use hits associated with current particle
          key_type trackId = hit.trackId();
          if( trackId != isp->first ) continue;

          // Get the hit information.

          unsigned int id = hit.volumeId();

          if( id<=0 || id>nvdet || ntp.isvd[id-1] ) continue;

          const CLHEP::Hep3Vector& pos = hit.position();
          const CLHEP::Hep3Vector& mom = hit.momentum();

          CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
          if( vdg->getRotation(id)!=0 ) {
            lpos *= *(vdg->getRotation(id));
          }

          ntp.isvd[id-1] = true;
          ntp.tvd[id-1]  = hit.time();
          ntp.gtvd[id-1] = hit.properTime();
          ntp.xvd[id-1]  = pos.x();
          ntp.yvd[id-1]  = pos.y();
          ntp.zvd[id-1]  = pos.z();
          ntp.pxvd[id-1] = mom.x();
          ntp.pyvd[id-1] = mom.y();
          ntp.pzvd[id-1] = mom.z();
          ntp.pvd[id-1]  = mom.mag();
          ntp.xlvd[id-1] = lpos.x();
          ntp.ylvd[id-1] = lpos.y();
          ntp.zlvd[id-1] = lpos.z();

        } // end loop over hits.

        // Keep only stopped particles
        if( _stopped_only && !ntp.isstop ) continue;

        // Keep only those particles which went through required VD
        if( _vd_required>0 && !ntp.isvd[_vd_required-1] ) continue;

        // Keep only those particles, which die late enough
        if( _timeCut>0.1 && ntp.tstop<_timeCut ) continue;

        _ntpart->Fill();

      }
    }

  }

}  // end namespace mu2e

using mu2e::ReadVirtualDetector;
DEFINE_ART_MODULE(ReadVirtualDetector);
