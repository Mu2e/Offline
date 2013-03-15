//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeCellsWithTurningOff_module.cc,v 1.4 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author G.F. Tassielli. Class derived by MakeStrawHit
//

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its tool chain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
//using art::Event;

namespace mu2e {

  // Utility class (structure) to hold calculated drift time for G4 hits

  class StepHit {
    typedef art::Handle<StepPointMCCollection> const* PHandle;
    
  public:
    
    art::Ptr<StepPointMC> _ptr;
    double _edep;
    double _dca;
    double _driftTime;
    double _distanceToMid;
    double _t1;
    
    StepHit( art::Ptr<StepPointMC> const& ptr, double edep, double dca, double driftT, double toMid, double t1):
      _ptr(ptr), _edep(edep), _dca(dca), _driftTime(driftT),
      _distanceToMid(toMid), _t1(t1) { }
    
    // This operator is overloaded in order to time-sort the hits
    bool operator <(const StepHit& b) const { return (_t1 < b._t1); }

    };

  //--------------------------------------------------------------------
  //
  //
  class MakeCellsWithTurningOff : public art::EDProducer {
  public:
    explicit MakeCellsWithTurningOff(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<string>("makerModuleLabel")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _poissonQ( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      // Control some information messages.
      _cellsToTurnOff(pset.get<int>("cellsToTurnOff",0)){

      // Tell the framework what we make.
      produces<StrawHitCollection>();
      produces<StrawHitMCTruthCollection>();
      produces<PtrStepPointMCVectorCollection>("StrawHitMCPtrWithTurningOff");

      _cellsToTurnOff = (_cellsToTurnOff<0?-_cellsToTurnOff: _poissonQ.fire(_cellsToTurnOff));


    }
    virtual ~MakeCellsWithTurningOff() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    typedef std::map<StrawIndex,std::vector<art::Ptr<StepPointMC> > > CellsWithTurningOffMap;

    // Diagnostics level.
    int _diagLevel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    std::string _makerModuleLabel;

    // Parameters
    std::string _g4ModuleLabel;  // Name of the module that made these hits.

    CLHEP::RandPoissonQ _poissonQ;

    int _cellsToTurnOff;

    int _nCells, _nSLayer;

    // A category for the error logger.
    const std::string _messageCategory;
    std::set< std::pair<int, int> > _offCellsIndex;

  };

  void MakeCellsWithTurningOff::beginJob(){

  }

  void
  MakeCellsWithTurningOff::produce(art::Event& event) {
    //cout<<"Event "<<event.event()<<endl;
    if ( _diagLevel > 0 ) cout << "MakeCellsWithTurningOff: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;
    art::ServiceHandle<GeometryService> geom;
    const Tracker& tracker = getTrackerOrThrow();
    const ITracker &itr = static_cast<const ITracker&>( tracker );
    CellGeometryHandle *itwp=itr.getCellGeometryHandle();

    if (ncalls == 1) {
      
      if (_diagLevel >0) cout << _cellsToTurnOff << " cells to be turned off\n" 
			      << "They are:" << endl;
      _nCells = itr.nSWire(); 
      _nSLayer = itr.nSuperLayers();
      
      art::RandomNumberGenerator::base_engine_t& engine = art::ServiceHandle<art::RandomNumberGenerator>()->getEngine();
      CLHEP::RandFlat flat(engine);;
      
      
      for (int i=0; i<_cellsToTurnOff; ++i) {
	int turnOffSLay = flat.fire() * _nSLayer;
	int turnOffCell = flat.fire() * _nCells;
	if (!(_offCellsIndex.insert( pair<int,int>(turnOffSLay, turnOffCell)).second)) {
	  i -= 1;
	} else {
	  if (_diagLevel>0) cout << "Superlayer " << turnOffSLay << "\tCell " << turnOffCell << endl;
	}
      }
      
      //check 
      if ((size_t)_cellsToTurnOff != _offCellsIndex.size()) {
	cout << "Something wrong in the routine to select the cells to be"
	     << " turned off" << endl;
      }
      
    }

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    event.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    // A container to hold the output hits.
    unique_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    unique_ptr<StrawHitMCTruthCollection>      truthHits(new StrawHitMCTruthCollection);
    unique_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);

    for (size_t i =0; i < hits->size(); ++i) {
      StrawHit const& hit(hits->at(i));
      StrawHitMCTruth const& MChit(hits_truth->at(i));
      PtrStepPointMCVector const& StepPointPtr(hits_mcptr->at(i));
      StrawIndex si = hit.strawIndex();
      //const Straw & str = tracker.getStraw(si);
      //const Cell & cell = static_cast<const Cell&>(str);
      itwp->SelectCellDet(si.asUint());
      int slayNo = itwp->GetSuperLayer();//cell.Id().getLayerId().getSuperLayer();
      int cellNo = itwp->GetWire();//cell.Id().getCell();
      if  (_offCellsIndex.find(pair<int,int>(slayNo,cellNo)) == _offCellsIndex.end()) {
	strawHits->push_back(StrawHit(hit));
	truthHits->push_back(StrawHitMCTruth(MChit));
	mcptrHits->push_back(PtrStepPointMCVector(StepPointPtr));
      }
    }
    event.put(std::move(strawHits));
    event.put(std::move(truthHits));
    event.put(std::move(mcptrHits),"StrawHitMCPtrWithTurningOff");
  }
}

using mu2e::MakeCellsWithTurningOff;
DEFINE_ART_MODULE(MakeCellsWithTurningOff);
