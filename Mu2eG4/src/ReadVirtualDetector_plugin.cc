//
// Plugin to read virtual detectors data and create ntuples
//
//  $Id: ReadVirtualDetector_plugin.cc,v 1.15 2011/05/17 15:36:00 greenc Exp $
//  $Author: greenc $
//  $Date: 2011/05/17 15:36:00 $
//
// Original author Ivan Logashenko
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadVirtualDetector : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;
    typedef SimParticleCollection::key_type key_type;

    explicit ReadVirtualDetector(fhicl::ParameterSet const& pset) : 
      _vdStepPoints(pset.get<string>("vdStepPoints","virtualdetector")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _ntvd(0), _ntpart(0) {
      
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

      nt = new float[200];

    }
  
    virtual ~ReadVirtualDetector() { }

    virtual void beginJob(art::EventSetup const&);
    virtual void beginRun(art::Run const&, art::EventSetup const& );

    void analyze(const art::Event& e, art::EventSetup const&);

  private:

    // Name of the VD StepPoint collection
    std::string _vdStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntvd;
    TTree* _ntpart;

    float *nt; // Need this buffer to fill TTree ntvd

    // Pointers to the physical volumes we are interested in
    // -- stopping target
    map<int,int> vid_stop;

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // List of virtual detectors to be saved
    set<int> vd_save;

  };
  
  void ReadVirtualDetector::beginJob(art::EventSetup const& ){

    vid_stop.clear();

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;
    
    _ntvd = tfs->make<TNtuple>( "ntvd", "Virtual Detectors ntuple", 
                                "evt:trk:sid:pdg:time:x:y:z:px:py:pz:xl:yl:zl:pxl:pyl:pzl:gtime");

    // Have to use TTree here, because one cannot use more than 100 variables in TNtuple

    _ntpart = tfs->make<TTree>("ntpart", "Particles ntuple");
    _ntpart->Branch("all",nt,
                    "evt:trk:pdg:"
                    "time:gtime:x:y:z:px:py:pz:"
                    "isstop:tstop:gtstop:xstop:ystop:zstop:"
                    "g4bl_evt:g4bl_trk:g4bl_weight:g4bl_time:"
                    "parent_id:parent_pdg:"
                    "parent_x:parent_y:parent_z:"
                    "parent_px:parent_py:parent_pz:"
                    "nvd:isvd[10]:"
                    "tvd[10]:gtvd[10]:xvd[10]:yvd[10]:zvd[10]:"
                    "pxvd[10]:pyvd[10]:pzvd[10]:"
                    "xlvd[10]:ylvd[10]:zlvd[10]"
                    );

  }

  void ReadVirtualDetector::beginRun(art::Run const& run, art::EventSetup const& ){

    // Get pointers to the physical volumes we are interested
    art::Handle<PhysicalVolumeInfoCollection> physVolumes;
    run.getByType(physVolumes);
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

  void ReadVirtualDetector::analyze(const art::Event& event, art::EventSetup const&) {

    ++_nAnalyzed;

    // Access virtual detectors geometry information
    // If not virtual detectors are defined, skip the rest

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel("g4run",_vdStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    art::Handle<G4BeamlineInfoCollection> g4beamlineData;
    event.getByType(g4beamlineData);
    bool haveG4BL = g4beamlineData.isValid();
    if ( haveG4BL ) haveG4BL = (g4beamlineData->size()==1);

    // Fill virtual detectors ntuple

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){
      
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

    // Fill tracks ntuple
    if( haveSimPart && pdg_save.size()>0 ) {

      const int nvdet = 10;
      const int id0 = 29;

      // Go through SimParticle container and analyze one particle at a time
      for ( SimParticleCollection::const_iterator isp=simParticles->begin();
            isp!=simParticles->end(); ++isp ){
        SimParticle const& sim = isp->second;

        // It particle PDG ID is not in the list - skip it
        if( pdg_save.find(sim.pdgId()) == pdg_save.end() ) continue;

        // Clean the buffer 
        for( int i=0; i<(id0+12*nvdet); ++i ) nt[i]=0;

        // Save SimParticle info
        nt[0] = event.id().event();    // event_id
        nt[1] = sim.id().asInt();      // track_id
        nt[2] = sim.pdgId();           // PDG id
        nt[3] = sim.startGlobalTime(); // start time
        nt[4] = sim.startProperTime(); // start time
        CLHEP::Hep3Vector const & pos_start = sim.startPosition();
        CLHEP::Hep3Vector const & mom_start = sim.startMomentum();
        nt[5] = pos_start.x();
        nt[6] = pos_start.y();
        nt[7] = pos_start.z();
        nt[8] = mom_start.x();
        nt[9] = mom_start.y();
        nt[10]= mom_start.z();

        // Check id of the volume there particle dies
        if( sim.endDefined() ) {
          if( vid_stop.find(sim.endVolumeIndex()) != vid_stop.end() ) nt[11] = 1;
          nt[12] = sim.endGlobalTime();
          nt[13] = sim.endProperTime();
          CLHEP::Hep3Vector const & pos_end = sim.endPosition();
          nt[14] = pos_end.x();
          nt[15] = pos_end.y();
          nt[16] = pos_end.z();
        }

        if( haveG4BL ) {
          G4BeamlineInfo const& extra = g4beamlineData->at(0);
          nt[17] = extra.eventId();
          nt[18] = extra.trackId();
          nt[19] = extra.weight();
          nt[20] = extra.time();
        }

        // Parent info
        if( sim.hasParent() ) {
          nt[21] = sim.parentId().asInt();
          SimParticle const* sim_parent = simParticles->findOrNull(sim.parentId());
          if( sim_parent ) {
            nt[22] = sim_parent->pdgId();
            CLHEP::Hep3Vector const & pos_parent = sim_parent->startPosition();
            CLHEP::Hep3Vector const & mom_parent = sim_parent->startMomentum();
            nt[23] = pos_parent.x();
            nt[24] = pos_parent.y();
            nt[25] = pos_parent.z();
            nt[26] = mom_parent.x();
            nt[27] = mom_parent.y();
            nt[28]= mom_parent.z();
          }
        } else {
          nt[21]=-1;
        }
          
        nt[id0] = nvdet;
      
        // Loop over all virtual detectors and fill corresponding data
        for ( size_t i=0; i<hits->size(); ++i ){

          // Alias, used for readability.
          const StepPointMC& hit = (*hits)[i];
          
          // Only use hits associated with current particle
          key_type trackId = hit.trackId();
          if( trackId != isp->first ) continue;

          // Get the hit information.
          
          int id = hit.volumeId();

          if( id<=0 || id>nvdet || (nt[id0+id]!=0) ) continue; 

          const CLHEP::Hep3Vector& pos = hit.position();
          const CLHEP::Hep3Vector& mom = hit.momentum();
 
          CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
          if( vdg->getRotation(id)!=0 ) {
            lpos *= *(vdg->getRotation(id));
          }

          nt[id0         +id] = 1.0;
          nt[id0+   nvdet+id] = hit.time();
          nt[id0+ 2*nvdet+id] = hit.properTime();
          nt[id0+ 3*nvdet+id] = pos.x();
          nt[id0+ 4*nvdet+id] = pos.y();
          nt[id0+ 5*nvdet+id] = pos.z();
          nt[id0+ 6*nvdet+id] = mom.x();
          nt[id0+ 7*nvdet+id] = mom.y();
          nt[id0+ 8*nvdet+id] = mom.z();
          nt[id0+ 9*nvdet+id] = lpos.x();
          nt[id0+10*nvdet+id] = lpos.y();
          nt[id0+11*nvdet+id] = lpos.z();

        } // end loop over hits.

        _ntpart->Fill();
        
      }
    }

  }

}  // end namespace mu2e

using mu2e::ReadVirtualDetector;
DEFINE_ART_MODULE(ReadVirtualDetector);
