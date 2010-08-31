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
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TH1F.h"
#include "TNtuple.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadVirtualDetector : public edm::EDAnalyzer {
  public:
    
    explicit ReadVirtualDetector::ReadVirtualDetector(edm::ParameterSet const& pset) : 
      _vdStepPoints(pset.getUntrackedParameter<string>("vdStepPoints","virtualdetector")),
      _ntup(0){
    }
  
    virtual ~ReadVirtualDetector() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Name of the VD StepPoint collection
    std::string _vdStepPoints;

    TNtuple* _ntup;

  };
  
  void ReadVirtualDetector::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;
    
    _ntup = tfs->make<TNtuple>( "ntvd", "Virtual Detectors ntuple", 
				"evt:trk:sid:pdg:time:x:y:z:px:py:pz:xl:yl:zl:pxl:pyl:pzl");

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

    // ntuple buffer.
    float nt[17];

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

      _ntup->Fill(nt);

    } // end loop over hits.

  }

}  // end namespace mu2e



using mu2e::ReadVirtualDetector;
DEFINE_FWK_MODULE(ReadVirtualDetector);
