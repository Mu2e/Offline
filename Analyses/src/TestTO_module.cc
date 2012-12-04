//
// A module to study background rates in the detector subsystems.
//
// $Id: TestTO_module.cc,v 1.2 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//
// Original author Gianni Onorato
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StatusG4.hh"
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {


  class TestTO : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;

    explicit TestTO(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _makerModuleLabelWTO(pset.get<std::string>("makerModuleLabelWTO")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run"))
    {
    }
    virtual ~TestTO() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    void doITracker(art::Event const& evt, bool skip);

    // Diagnostic level
    int _diagLevel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;
    std::string _makerModuleLabelWTO;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

  };

  void TestTO::beginJob( ) {
  }

  void TestTO::analyze(art::Event const& evt ) {

    doITracker(evt, false);
    
  } // end of analyze

  void TestTO::endJob() {
  }
    
  void TestTO::doITracker(art::Event const& evt, bool skip) {

    if (skip) return;

    const Tracker& tracker = getTrackerOrThrow();

    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    if (!(hits->size() == hits_truth->size() &&
          hits_mcptr->size() == hits->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits: " << hits->size()
        << " MCTruthStrawHits: " << hits_truth->size()
        << " MCPtr: " << hits_mcptr->size();
    }

    art::Handle<StrawHitCollection> pdataHandleWTO;
    evt.getByLabel(_makerModuleLabelWTO,pdataHandleWTO);
    StrawHitCollection const* hitsWTO = pdataHandleWTO.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandleWTO;
    evt.getByLabel(_makerModuleLabelWTO,truthHandleWTO);
    StrawHitMCTruthCollection const* hits_truthWTO = truthHandleWTO.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandleWTO;
    evt.getByLabel(_makerModuleLabelWTO,"StrawHitMCPtrWithTurningOff",mcptrHandleWTO);
    PtrStepPointMCVectorCollection const* hits_mcptrWTO = mcptrHandleWTO.product();

    if (!(hitsWTO->size() == hits_truthWTO->size() &&
          hits_mcptrWTO->size() == hitsWTO->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits with Turning Off: " << hitsWTO->size()
        << " MCTruthStrawHits with Turning Off: " << hits_truthWTO->size()
        << " MCPtr with Turning Off: " << hits_mcptrWTO->size();
    }

    size_t nStrawPerEvent = hits->size();
    size_t nStrawPerEventWTO = hitsWTO->size();

    if (nStrawPerEvent != nStrawPerEventWTO) {
      cout << "************* Something has been turned off ************\n\n\n" << endl;
    }
    
    for (size_t i=0; i<nStrawPerEvent; ++i) {

      StrawHit             const&  hit(hits->at(i));
      StrawHitMCTruth      const&  truth(hits_truth->at(i));
      PtrStepPointMCVector const&  mcptr(hits_mcptr->at(i));

      StrawIndex si = hit.strawIndex();
      const Straw & str = tracker.getStraw(si);
      const Cell & cell = static_cast<const Cell&>( str );

      int sid = cell.Id().getCell();
      int lid = cell.Id().getLayer();
      int did = cell.Id().getLayerId().getSuperLayer();
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      cout << "Cell: " << sid << " , " << lid << " , " << did << " : " << endl;
      cout << "\t\tdriftTime = " << driftTime << " - driftDistance = "
	   << driftDistance << endl;
      cout << "\t\tThere are " << mcptr.size() << " stepPoints..." << endl; 

      for (size_t j = 0; j < mcptr.size(); ++j) {
	StepPointMC const& mchit = *mcptr[j];
	double hitE = mchit.eDep();
	cout << "\t\t\tStepPoint " << j+1 << " - eDep = " << hitE << endl;
      }
      
    } //end of Cellhits loop


    for (size_t i=0; i<nStrawPerEventWTO; ++i) {

      StrawHit             const&  hit(hitsWTO->at(i));
      StrawHitMCTruth      const&  truth(hits_truthWTO->at(i));
      PtrStepPointMCVector const&  mcptr(hits_mcptrWTO->at(i));

      StrawIndex si = hit.strawIndex();
      const Straw & str = tracker.getStraw(si);
      const Cell & cell = static_cast<const Cell&>( str );

      int sid = cell.Id().getCell();
      int lid = cell.Id().getLayer();
      int did = cell.Id().getLayerId().getSuperLayer();
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      cout << "Cell with turning off: " << sid << " , " << lid << " , " << did << " : " << endl;
      cout << "\t\tdriftTime = " << driftTime << " - driftDistance = "
	   << driftDistance << endl;
      cout << "\t\tThere are " << mcptr.size() << " stepPoints..." << endl; 

      for (size_t j = 0; j < mcptr.size(); ++j) {
	StepPointMC const& mchit = *mcptr[j];
	double hitE = mchit.eDep();
	cout << "\t\t\tStepPoint " << j+1 << " - eDep = " << hitE << endl;
      }
      
    } //end of Cellhits loop

  } // end of doITracker
  
}

using mu2e::TestTO;
DEFINE_ART_MODULE(TestTO);

