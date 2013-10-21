//
// An EDAnalyzer module that reads back ExtMonUCI hits created by G4 and makes trees.
//
// $Id: ExtMonUCIReadBack_module.cc,v 1.3 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "RecoDataProducts/inc/ExtMonUCITofHitCollection.hh"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
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

  class ExtMonUCIReadBack : public art::EDAnalyzer {
  public:

    explicit ExtMonUCIReadBack(fhicl::ParameterSet const& pset);
    virtual ~ExtMonUCIReadBack() { }

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);

  private:

    // Start: run time parameters

    // Diagnostics printout level
    int _diagLevel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Module label of the generator module that was passed as input to G4.
    std::string _generatorModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Module which made the ExtMonUCIHits
    std::string _extMonUCIModuleLabel;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Limit the size of the TGraph.
    int _xyHitsMax;

    // End: run time parameters

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms, ntuples, TGraphs.
    // ExtMonUCI ntuples
    TNtuple* _ntupExtMonUCITof;
    TNtuple* _ntupExtMonUCITofMC;

    int _nBadG4Status;

    // Do the work specific to one of the trackers.
    void doExtMonUCI(const art::Event& event);

  };

  ExtMonUCIReadBack::ExtMonUCIReadBack(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _extMonUCIModuleLabel(pset.get<string>("extMonUCITofModuleLabel","ExtMonUCITofHitsMaker")),
    _minimumEnergy(pset.get<double>("minimumEnergy")),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _xyHitsMax(pset.get<int>("xyHitsMax",10000)),

    // Histograms
    _nAnalyzed(0),
    _ntupExtMonUCITof(0),
    _ntupExtMonUCITofMC(0),
    // Remaining member data
    _nBadG4Status(0){
  }

  void ExtMonUCIReadBack::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // ExtMonUCI Tof ntuple.
    _ntupExtMonUCITof = tfs->make<TNtuple>( "ntupExtMonUCITof", "Extinction Monitor UCI Tof Hits",
                                   "run:evt:stId:segId:t:edep:x:y:z");
    
    _ntupExtMonUCITofMC = tfs->make<TNtuple>( "ntupExtMonUCITofMC", "Extinction Monitor UCI Tof Hits MC Truth",
                                   "run:evt:stId:segId:t:edep:trk:pdgId:x:y:z:px:py:pz:vx:vy:vz:vpx:vpy:vpz:vt:primary:otrk:opdgId:ox:oy:oz:opx:opy:opz:ot");
    
  }

  void ExtMonUCIReadBack::analyze(const art::Event& event) {

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Inquire about the completion status of G4.
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    if ( _nAnalyzed < _maxFullPrint ){
      cerr << g4Status << endl;
    }

    // Abort if G4 did not complete correctly.
    // Use your own judgement about whether to abort or to continue.
    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting ExtMonUCIReadBack::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    // Call code appropriate for the tracker that is installed in this job.
    art::ServiceHandle<GeometryService> geom;

    if(geom->hasElement<ExtMonUCI::ExtMon>() ) {
      doExtMonUCI(event);
    }

  }

  void ExtMonUCIReadBack::doExtMonUCI(const art::Event& event) {

    // Gometry for the extMonUCI.
    GeomHandle<ExtMonUCI::ExtMon> extMonUCI;
    int nTofStations = extMonUCI->nTofStations();
    int nTofSegments = extMonUCI->nTofSegments();
    if ( _diagLevel > 1 && _nAnalyzed < _maxFullPrint ){
      cout << "Readback: ExtMonUCI nTofStations " << nTofStations << " nTofSegments " << nTofSegments << endl;
    }

    // Get handles to extmonuci collections
    art::Handle<ExtMonUCITofHitCollection> tofHits;
    art::Handle<ExtMonUCITofHitMCTruthCollection> tofMC;
    event.getByLabel(_extMonUCIModuleLabel,tofHits);
    event.getByLabel(_extMonUCIModuleLabel,tofMC);
          
    // Find pointers to the original G4 steps
    art::Handle<PtrStepPointMCVectorCollection> tofMCPtr;
    event.getByLabel(_extMonUCIModuleLabel,tofMCPtr);
      
    bool haveExtMonUCITof = ( tofHits.isValid() && tofMC.isValid() );
      
    if( ! haveExtMonUCITof) return;
      
    //art::ServiceHandle<GeometryService> geom;
    //if( ! geom->hasElement<Calorimeter>() ) return;
    //GeomHandle<Calorimeter> cg;

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){    
      cout << "EVENT " << event.id().run() << " subRunID " << event.id().subRun() << " event " << event.id().event() << endl; 
      cout << " time " << event.time().value() << endl;
    }

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<tofHits->size(); ++i ) {
        ExtMonUCITofHit const & hit = (*tofHits).at(i);
        cout << "Readback: ExtMonUCITofHit " << hit << endl;
      }
    } 
      
    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<tofMC->size(); ++i ) {
        ExtMonUCITofHitMCTruth const & hit = (*tofMC).at(i);
        cout << "Readback: ExtMonUCITofHitMCTruth " << hit << endl;
      }
    } 
      
      
    // ntuple buffer.
    float nt[_ntupExtMonUCITof->GetNvar()];
      
    // Loop over all hits.
    for ( size_t i=0; i<tofHits->size(); ++i ){

      // Alias, used for readability.
      const ExtMonUCITofHit& hit = (*tofHits)[i];

      // Fill the ntuple.
      nt[ 0] = event.id().run();
      nt[ 1] = event.id().event();
      nt[ 2] = hit.stationId();
      nt[ 3] = hit.segmentId();    
      nt[ 4] = hit.time();
      nt[ 5] = hit.energyDep();       
      CLHEP::Hep3Vector const &  mid = extMonUCI->tof(hit.stationId(), hit.segmentId())->origin();
      nt[ 6] = mid.x();
      nt[ 7] = mid.y();
      nt[ 8] = mid.z();
             
      _ntupExtMonUCITof->Fill(nt); 
             
    }        
            
    // ntuple buffer.              
    float ntmc[_ntupExtMonUCITofMC->GetNvar()];
             
    // Loop over all hits.
    for ( size_t i=0; i<tofMC->size(); ++i ){

      // Alias, used for readability.
      const ExtMonUCITofHitMCTruth& hit = (*tofMC)[i];
  
      // Fill the ntuple.
      ntmc[ 0] = event.id().run();
      ntmc[ 1] = event.id().event();
      ntmc[ 2] = hit.stationId();
      ntmc[ 3] = hit.segmentId();
      ntmc[ 4] = hit.time();
      ntmc[ 5] = hit.energyDep();
      ntmc[ 6] = hit.trackId();
      ntmc[ 7] = hit.pdgId();
      ntmc[ 8] = hit.position().x();
      ntmc[ 9] = hit.position().y();
      ntmc[10] = hit.position().z();
      ntmc[11] = hit.momentum().x();
      ntmc[12] = hit.momentum().y();
      ntmc[13] = hit.momentum().z();
      ntmc[14] = hit.vertex().x();
      ntmc[15] = hit.vertex().y();
      ntmc[16] = hit.vertex().z();
      ntmc[17] = hit.vertexMomentum().x();
      ntmc[18] = hit.vertexMomentum().y();
      ntmc[19] = hit.vertexMomentum().z();
      ntmc[20] = hit.vertexTime();
      ntmc[21] = hit.isPrimary();
      ntmc[22] = hit.orgTrackId();
      ntmc[23] = hit.orgPdgId();
      ntmc[24] = hit.orgVertex().x();
      ntmc[25] = hit.orgVertex().y();
      ntmc[26] = hit.orgVertex().z();
      ntmc[27] = hit.orgVertexMomentum().x();
      ntmc[28] = hit.orgVertexMomentum().y();
      ntmc[29] = hit.orgVertexMomentum().z();
      ntmc[30] = hit.orgTime();

      _ntupExtMonUCITofMC->Fill(ntmc);

    }

  } // end of doExtMonUCI

  void ExtMonUCIReadBack::endJob(){
    cout << "ExtMonUCIReadBack::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
         << endl;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonUCIReadBack);
