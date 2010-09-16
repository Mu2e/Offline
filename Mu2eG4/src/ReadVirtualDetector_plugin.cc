//
// Plugin to read virtual detectors data and create ntuples
//
// Original author Ivan Logashenko
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
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

  class ReadVirtualDetector : public edm::EDAnalyzer {
  public:

    typedef vector<int> vint;

    explicit ReadVirtualDetector::ReadVirtualDetector(edm::ParameterSet const& pset) : 
      _vdStepPoints(pset.getUntrackedParameter<string>("vdStepPoints","virtualdetector")),
      _ntvd(0), _ntpart(0) {
      
      vint const & pdg_ids = pset.getUntrackedParameter<vint>("savePDG", vint());
      if( pdg_ids.size()>0 ) {
	cout << "ReadVirtualDetector: save following particle types in the ntuple: ";
	for( int i=0; i<pdg_ids.size(); ++i ) {
	  pdg_save.insert(pdg_ids[i]);
	  cout << pdg_ids[i] << ", ";
	}
	cout << endl;
      }

      nt = new float[200];

    }
  
    virtual ~ReadVirtualDetector() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual void beginRun(edm::Run const&, edm::EventSetup const& );

    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Name of the VD StepPoint collection
    std::string _vdStepPoints;

    TNtuple* _ntvd;
    TTree* _ntpart;

    float *nt; // Need this buffer to fill TTree
    
    // Pointers to the physical volumes we are interested in
    // -- stopping target
    map<int,int> vid_stop;

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

  };
  
  void ReadVirtualDetector::beginJob(edm::EventSetup const& ){

    vid_stop.clear();

    // Get access to the TFile service.

    edm::Service<edm::TFileService> tfs;
    
    _ntvd = tfs->make<TNtuple>( "ntvd", "Virtual Detectors ntuple", 
				"evt:trk:sid:pdg:time:x:y:z:px:py:pz:xl:yl:zl:pxl:pyl:pzl:gtime");

    // Have to use TTree here, because one cannot use more than 100 variables in TNtuple

    _ntpart = tfs->make<TTree>("ntpart", "Particles ntuple");
    _ntpart->Branch("all",nt,
		    "evt:trk:pdg:"
		    "time:gtime:x:y:z:px:py:pz:"
		    "isstop:tstop:gtstop:xstop:ystop:zstop:"
		    "g4bl_evt:g4bl_trk:g4bl_weight:g4bl_time:"
		    "nvd:isvd[10]:"
		    "tvd[10]:gtvd[10]:xvd[10]:yvd[10]:zvd[10]:"
		    "pxvd[10]:pyvd[10]:pzvd[10]:"
		    "xlvd[10]:ylvd[10]:zlvd[10]");

  }

  void ReadVirtualDetector::beginRun(edm::Run const& run, edm::EventSetup const& ){

    // Get pointers to the physical volumes we are interested
    edm::Handle<PhysicalVolumeInfoCollection> physVolumes;
    run.getByType(physVolumes);
    if( physVolumes.isValid() ) {

      for ( int i=0; i<physVolumes->size(); ++i ) {
	if( (*physVolumes)[i].name() == "TargetFoil_" ) {
	  vid_stop[i] = (*physVolumes)[i].copyNo();
	  cout << "ReadVirtualDetector: register stopping target volume " << i << " = "
	       << (*physVolumes)[i].name() << " " << (*physVolumes)[i].copyNo() << endl;
	}
      }

    }

  }

  void ReadVirtualDetector::analyze(const edm::Event& event, edm::EventSetup const&) {

    // Access virtual detectors geometry information
    // If not virtual detectors are defined, skip the rest

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel("g4run",_vdStepPoints,hits);

    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    edm::Handle<G4BeamlineInfoCollection> g4beamlineData;
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

      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();
 
      CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
      CLHEP::Hep3Vector lmom = mom;
      if( vdg->getRotation(id)!=0 ) {
	lpos *= *(vdg->getRotation(id));
	lmom *= *(vdg->getRotation(id));
      }

      // Get track info
      int trackId = hit.trackId();
      int pdgId = 0;
      if ( haveSimPart ){
        SimParticle const& sim = simParticles->at(trackId);
        pdgId = sim.pdgId();
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId;
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

    } // end loop over hits.

    // Fill tracks ntuple
    if( haveSimPart && pdg_save.size()>0 ) {

      const int nvdet = 10;
      const int id0 = 21;

      // Go through SimParticle container and analyze one particle at a time
      for( int isp=0; isp<simParticles->size(); ++isp ) {

        SimParticle const& sim = simParticles->at(isp);

	// It particle PDG ID is not in the list - skip it
	if( pdg_save.find(sim.pdgId()) == pdg_save.end() ) continue;
	bool toSave = false;

	// Clean the buffer 
	for( int i=0; i<(id0+12*nvdet); ++i ) nt[i]=0;

	// Save SimParticle info
	nt[0] = event.id().event();    // event_id
	nt[1] = isp;                   //track_id
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

	nt[id0] = nvdet;
      
	// Loop over all virtual detectors and fill corresponding data
	for ( size_t i=0; i<hits->size(); ++i ){

	  // Alias, used for readability.
	  const StepPointMC& hit = (*hits)[i];
	  
	  // Only use hits associated with current particle
	  int trackId = hit.trackId();
	  if( trackId != isp ) continue;

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
DEFINE_FWK_MODULE(ReadVirtualDetector);
