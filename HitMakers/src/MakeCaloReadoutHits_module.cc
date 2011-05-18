//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalOnlyHit, CrystalHitMCPtr
// objects.
//
// Original author Ivan Logashenko.
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalOnlyHitCollection.hh"
#include "ToyDP/inc/DPIndexVector.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/sort_functors.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using art::Event;

namespace mu2e {

  // Utility class (structure) to hold RO info for G4 hits

  class ROHit {
  public:

    int    _hit_id;
    double _edep;
    double _edep_corr;
    int    _charged;
    double _time;

    ROHit(int hit_id, double edep, double edep1, int charged, double time):
      _hit_id(hit_id), _edep(edep), _edep_corr(edep1),
      _charged(charged), _time(time) { }

    // This operator is overloaded in order to time-sort the hits
    bool operator <(const ROHit& b) const { return (_time < b._time); }

  };

  //--------------------------------------------------------------------
  //
  //
  class MakeCaloReadoutHits : public art::EDProducer {
  public:
    explicit MakeCaloReadoutHits(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _stepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
      _rostepPoints(pset.get<string>("calorimeterROStepPoints","calorimeterRO")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),

      _messageCategory("CaloReadoutHitsMaker"){

      // Tell the framework what we make.
      produces<CaloHitCollection>();
      produces<CaloHitMCTruthCollection>();
      produces<CaloCrystalOnlyHitCollection>();
      produces<DPIndexVectorCollection>("CaloHitMCCrystalPtr");
      produces<DPIndexVectorCollection>("CaloHitMCReadoutPtr");

    }
    virtual ~MakeCaloReadoutHits() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StepPoint collection
    std::string _stepPoints;
    std::string _rostepPoints;

    // Parameters
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // A category for the error logger.
    const std::string _messageCategory;

    void makeCalorimeterHits (const art::Handle<StepPointMCCollection>&,
			      const art::Handle<StepPointMCCollection>&,
			      CaloHitCollection &,
			      CaloHitMCTruthCollection&,
			      CaloCrystalOnlyHitCollection&,
			      DPIndexVectorCollection&,
			      DPIndexVectorCollection&);

  };

  void MakeCaloReadoutHits::beginJob(){

  }

  void
  MakeCaloReadoutHits::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;

    // A container to hold the output hits.
    auto_ptr<CaloHitCollection>               caloHits(new CaloHitCollection);
    auto_ptr<CaloHitMCTruthCollection>        caloMCHits(new CaloHitMCTruthCollection);
    auto_ptr<CaloCrystalOnlyHitCollection> caloCrystalMCHits(new CaloCrystalOnlyHitCollection);
    auto_ptr<DPIndexVectorCollection>         caloMCptrHits(new DPIndexVectorCollection);
    auto_ptr<DPIndexVectorCollection>         caloMCroptrHits(new DPIndexVectorCollection);

    // Ask the event to give us a handle to the requested hits.

    art::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_stepPoints,points);
    int nHits = points->size();

    art::Handle<StepPointMCCollection> rohits;
    event.getByLabel(_g4ModuleLabel,_rostepPoints,rohits);
    int nroHits = rohits->size();

    if( nHits>0 || nroHits>0 ) {
      makeCalorimeterHits(points, rohits,
			  *caloHits, *caloMCHits, *caloCrystalMCHits,
			  *caloMCptrHits, *caloMCroptrHits);
    }

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeCaloReadoutHits: Total number of calorimeter hits = "
	   << caloHits->size()
	   << endl;
      cout << "MakeCaloReadoutHits: Total number of crystal MC hits = "
	   << caloCrystalMCHits->size()
	   << endl;
    }

    // Add the output hit collection to the event
    event.put(caloHits);
    event.put(caloMCHits);
    event.put(caloCrystalMCHits);
    event.put(caloMCptrHits,"CaloHitMCCrystalPtr");
    event.put(caloMCroptrHits,"CaloHitMCReadoutPtr");

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() end" << endl;

  } // end of ::analyze.

  void MakeCaloReadoutHits::makeCalorimeterHits (const art::Handle<StepPointMCCollection>& steps,
				      const art::Handle<StepPointMCCollection>& rosteps,
				      CaloHitCollection &caloHits,
				      CaloHitMCTruthCollection& caloHitsMCTruth,
				      CaloCrystalOnlyHitCollection& caloCrystalHitsMCTruth,
				      DPIndexVectorCollection& caloHitsMCCrystalPtr,
				      DPIndexVectorCollection& caloHitsMCReadoutPtr ) {

    // Product Id of the input points.
    art::ProductID const& crystal_id(steps.id());
    art::ProductID const& readout_id(rosteps.id());

    // Get calorimeter geometry description
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    double length     = cal.crystalHalfLength();
    double nonUniform = cal.getNonuniformity();
    double timeGap    = cal.getTimeGap();
    double addEdep    = cal.getElectronEdep();
    int    nro        = cal.nROPerCrystal();

    // Organize steps by readout elements

    int nSteps   = steps->size();
    int nroSteps = rosteps->size();

    // First vector is list of crystal steps, associated with parcular readout element.
    // Second vector is list of readout steps, associated with parcular readout element.
    typedef std::map<int,std::pair<std::vector<int>,std::vector<int> > > HitMap;
    HitMap hitmap;

    for ( int i=0; i<nSteps; ++i){
      StepPointMC const& h = (*steps)[i];
      for( int j=0; j<nro; ++j ) {
	vector<int> &steps_id = hitmap[h.volumeId()*nro+j].first;
	steps_id.push_back(i);
      }
    }

    for ( int i=0; i<nroSteps; ++i){
      StepPointMC const& h = (*rosteps)[i];
      vector<int> &steps_id = hitmap[h.volumeId()].second;
      steps_id.push_back(i);
    }

    // Loop over all readout elements to form ro hits

    vector<ROHit> ro_hits;

    for(HitMap::const_iterator ro = hitmap.begin(); ro != hitmap.end(); ++ro ) {

      // readout ID
      int roid = ro->first;

      // Prepare info for hit creation before the next iteration
      ro_hits.clear();

      // Loop over all hits found for this readout element

      vector<int> const& isteps = ro->second.first;
      vector<int> const& irosteps = ro->second.second;

      // Loop over steps inside the crystal

      for( size_t i=0; i<isteps.size(); i++ ) {

        int hitRef = isteps[i];
	StepPointMC const& h = (*steps)[hitRef];

        double edep    = h.eDep()/nro; // each ro has its crystal energy deposit assigned to it
	if( edep<=0.0 ) continue; // Do not create hit if there is no energy deposition

	// Hit position in Mu2e frame
        CLHEP::Hep3Vector const& pos = h.position();
	// Hit position in local crystal frame
        CLHEP::Hep3Vector posLocal = cal.toCrystalFrame(roid,pos);
	// Calculate correction for edep
	double edep_corr = edep * (1.0+(posLocal.z()/length)*nonUniform/2.0);

	ro_hits.push_back(ROHit(hitRef,edep,edep_corr,0,h.time()));

      }

      // Loop over steps inside the readout (direct energy deposition in APD)

      for( size_t i=0; i<irosteps.size(); i++ ) {

        int hitRef = irosteps[i];
	StepPointMC const& h = (*rosteps)[hitRef];

	// There is no cut on energy deposition here - may be, we need to add one?

        double hitTime = h.time();

	ro_hits.push_back(ROHit(hitRef,0,0,1,hitTime));

      }

      // Sort hits by time
      sort(ro_hits.begin(), ro_hits.end());

      // Loop over sorted hits and form complete ro/calorimeter hits

      double h_time    = ro_hits[0]._time;
      double h_edep    = ro_hits[0]._edep;
      double h_edepc   = ro_hits[0]._edep_corr;
      int    h_charged = ro_hits[0]._charged;
      DPIndexVector mcptr_crystal;
      DPIndexVector mcptr_readout;
      if( ro_hits[0]._charged==0 ) {
	mcptr_crystal.push_back(DPIndex(crystal_id,ro_hits[0]._hit_id));
      } else {
	mcptr_readout.push_back(DPIndex(readout_id,ro_hits[0]._hit_id));
      }

      for( size_t i=1; i<ro_hits.size(); ++i ) {
	if( (ro_hits[i]._time-ro_hits[i-1]._time) > timeGap ) {
	  // Save current hit
	  caloHits.push_back(       CaloHit(       roid,h_time,h_edepc+h_charged*addEdep));
	  caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_charged));
	  caloHitsMCCrystalPtr.push_back(mcptr_crystal);
	  caloHitsMCReadoutPtr.push_back(mcptr_readout);
	  // ...and create new hit
	  mcptr_crystal.clear();
	  mcptr_readout.clear();
	  if( ro_hits[i]._charged==0 ) {
	    mcptr_crystal.push_back(DPIndex(crystal_id,ro_hits[i]._hit_id));
	  } else {
	    mcptr_readout.push_back(DPIndex(readout_id,ro_hits[i]._hit_id));
	  }
	  h_time    = ro_hits[i]._time;
	  h_edep    = ro_hits[i]._edep;
	  h_edepc   = ro_hits[i]._edep_corr;
	  h_charged = ro_hits[i]._charged;
	} else {
	  // Append data to hit
	  h_edep  += ro_hits[i]._edep;
	  h_edepc += ro_hits[i]._edep_corr;
	  if( ro_hits[i]._charged>0 ) h_charged = 1; // this does not count the charge...
	  if( ro_hits[i]._charged==0 ) {
	    mcptr_crystal.push_back(DPIndex(crystal_id,ro_hits[i]._hit_id));
	  } else {
	    mcptr_readout.push_back(DPIndex(readout_id,ro_hits[i]._hit_id));
	  }
	}
      }

      /*
      cout << "CaloHit: id=" << roid
	   << " time=" << h_time
	   << " trueEdep=" << h_edep
	   << " corrEdep=" << h_edepc
	   << " chargedEdep=" << (h_edepc+h_charged*addEdep)
	   << endl;
      */

      caloHits.push_back(       CaloHit(roid,h_time,h_edepc+h_charged*addEdep));
      caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_charged));
      caloHitsMCCrystalPtr.push_back(mcptr_crystal);
      caloHitsMCReadoutPtr.push_back(mcptr_readout);

    }

    // now CaloCrystalOnlyHit
    // Organize hits by crystal id; reject all but first RO; reject all RO hit directly

    hitmap.clear();

    for ( int i=0; i<nSteps; ++i){
      StepPointMC const& h = (*steps)[i];
      hitmap[h.volumeId()].first.push_back(i);
    }

    CaloCrystalOnlyHitCollection cr_hits;

    for(HitMap::const_iterator cr = hitmap.begin(); cr != hitmap.end(); ++cr) {

      // crystal id
      int cid = cr->first;

      // Prepare info for hit creation before the next iteration
      cr_hits.clear();

      // Loop over all hits found for this crystal
      vector<int> const& isteps = cr->second.first;

      for( size_t i=0; i<isteps.size(); i++ ) {
	StepPointMC const& h = (*steps)[isteps[i]];
	cr_hits.push_back(CaloCrystalOnlyHit(cid,h.time(),h.eDep()));
      }

      sort(cr_hits.begin(), cr_hits.end(), lessByTime<CaloCrystalOnlyHitCollection::value_type>());

      // now form final hits if they are close enough in time

      CaloCrystalOnlyHitCollection::value_type cHitMCTruth  = cr_hits[0];

      for ( size_t i=1; i<cr_hits.size(); ++i ) {

	if ( (cr_hits[i].time()-cr_hits[i-1].time()) > timeGap ) {

	  // Save current hit

          caloCrystalHitsMCTruth.push_back(cHitMCTruth);

	  // ...and create new hit

          cHitMCTruth  = cr_hits[i];

	} else {

	  // Add energy to the old hit (keep the "earlier" time)

          cHitMCTruth.setEnergyDep(cHitMCTruth.energyDep()+cr_hits[i].energyDep());

	}

      }

      caloCrystalHitsMCTruth.push_back(cHitMCTruth);

    }


  }

}

using mu2e::MakeCaloReadoutHits;
DEFINE_ART_MODULE(MakeCaloReadoutHits);
