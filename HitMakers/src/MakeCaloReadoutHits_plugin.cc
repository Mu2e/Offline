//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalHitMCTruth, CrystalHitMCPtr
// objects.
//
// Original author Ivan Logashenko.
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalHitMCTruthCollection.hh"
#include "Mu2eUtilities/inc/sort_functors.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using edm::Event;

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
    bool ROHit::operator <(const ROHit& b) const { return (_time < b._time); }

  };

  //--------------------------------------------------------------------
  //
  // 
  class MakeCaloReadoutHits : public edm::EDProducer {
  public:
    explicit MakeCaloReadoutHits(edm::ParameterSet const& pset) : 

      // Parameters
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _stepPoints(pset.getUntrackedParameter<string>("calorimeterStepPoints","calorimeter")),
      _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),

      _messageCategory("CaloReadoutHitsMaker"){

      // Tell the framework what we make.
      produces<CaloHitCollection>();
      produces<CaloHitMCTruthCollection>();
      produces<CaloCrystalHitMCTruthCollection>();

    }
    virtual ~MakeCaloReadoutHits() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void produce( edm::Event& e, edm::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StepPoint collection
    std::string _stepPoints;

    // Parameters
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // A category for the error logger.
    const std::string _messageCategory;

    void makeCalorimeterHits (const edm::Handle<StepPointMCCollection>&, 
			      CaloHitCollection &, 
			      CaloHitMCTruthCollection&,
			      CaloCrystalHitMCTruthCollection&);

  };

  void MakeCaloReadoutHits::beginJob(edm::EventSetup const& ){
    
  }

  void
  MakeCaloReadoutHits::produce(edm::Event& event, edm::EventSetup const&) {

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() begin" << endl;
      
    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    edm::Service<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;

    // A container to hold the output hits.
    auto_ptr<CaloHitCollection>               caloHits(new CaloHitCollection);
    auto_ptr<CaloHitMCTruthCollection>        caloMCHits(new CaloHitMCTruthCollection);
    auto_ptr<CaloCrystalHitMCTruthCollection> caloCrystalMCHits(new CaloCrystalHitMCTruthCollection);

    // Ask the event to give us a handle to the requested hits.
    edm::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_stepPoints,points);
    int nHits = points->size();

    // Product Id of the input points.
    edm::ProductID const& id(points.id());

    if( nHits>0 ) {
      makeCalorimeterHits(points, *caloHits, *caloMCHits, *caloCrystalMCHits);
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

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() end" << endl;

  } // end of ::analyze.
 
  void MakeCaloReadoutHits::makeCalorimeterHits (const edm::Handle<StepPointMCCollection>& hits, 
			     CaloHitCollection &caloHits, 
			     CaloHitMCTruthCollection& caloHitsMCTruth,
			     CaloCrystalHitMCTruthCollection& caloCrystalHitsMCTruth) {

    // Get calorimeter geometry description
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
    
    double length     = cal.crystalHalfLength();
    double nonUniform = cal.getNonuniformity();
    double timeGap    = cal.getTimeGap();
    double addEdep    = cal.getElectronEdep();
    int    nro        = cal.nROPerCrystal();

    // Organize hits by readout elements

    int nHits = hits->size();

    typedef std::map<int,std::vector<int> > HitMap;
    HitMap hitmap;
    for ( int i=0; i<nHits; ++i){
      StepPointMC const& h = (*hits)[i];
      vector<int> &hits_id = hitmap[h.volumeId()];
      hits_id.push_back(i);
    }

    // Loop over all readout elements to form ro hits

    vector<ROHit> ro_hits;

    for(HitMap::const_iterator ro = hitmap.begin(); ro != hitmap.end(); ++ro ) {

      // readout ID
      int roid = ro->first;

      // Prepare info for hit creation before the next iteration
      ro_hits.clear();

      // Loop over all hits found for this readout element

      vector<int> const& ihits = ro->second;

      for( size_t i=0; i<ihits.size(); i++ ) {

        int hitRef = ihits[i];
	StepPointMC const& h = (*hits)[hitRef];
        CLHEP::Hep3Vector const& pos = h.position();
        double edep    = h.eDep()/nro; // each ro has its crystal energy deposit assigned to it
        double hitTime = h.time();
	int charged = 0;

	// Calculate correction for edep; edep<0 means energy deposition in APD, not in crystal

	double edep_corr=0;
	if( edep>0.0 ) {
	  edep_corr = edep * (1.0+(pos.z()/length)*nonUniform/2.0);
        } else if ( edep==0.0 ) {
          if ( pos == CLHEP::Hep3Vector(0.0,0.0,0.0) ) {
            // it is most likely directly hit ro with no energy deposit, should it be "charged" ?
            charged=1;
          }
	} else { 
	  charged=1;
	  edep=0.0; 
	}
	
	ro_hits.push_back(ROHit(hitRef,edep,edep_corr,charged,hitTime));

      }

      // Sort hits by time
      sort(ro_hits.begin(), ro_hits.end());

      // Loop over sorted hits and form complete ro/calorimeter hits

      double h_time    = ro_hits[0]._time;
      double h_edep    = ro_hits[0]._edep;
      double h_edepc   = ro_hits[0]._edep_corr;
      int    h_charged = ro_hits[0]._charged;

      for( size_t i=1; i<ro_hits.size(); ++i ) {
	if( (ro_hits[i]._time-ro_hits[i-1]._time) > timeGap ) {
	  // Save current hit
	  caloHits.push_back(       CaloHit(       roid,h_time,h_edepc+h_charged*addEdep));
	  caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_charged));
	  // ...and create new hit	  
	  h_time    = ro_hits[i]._time;
	  h_edep    = ro_hits[i]._edep;
	  h_edepc   = ro_hits[i]._edep_corr;
	  h_charged = ro_hits[i]._charged;
	} else {
	  // Append data to hit
	  h_edep  += ro_hits[i]._edep;
	  h_edepc += ro_hits[i]._edep_corr;
	  if( ro_hits[i]._charged>0 ) h_charged = 1; // this does not count the charge...
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
      
    }

    // now CaloCrystalHitMCTruth
    // Organize hits by crystal id; reject all but first RO; reject all RO hit directly

    hitmap.clear();

    for ( int i=0; i<nHits; ++i){
      StepPointMC const& h = (*hits)[i];
      // reject all "non" first RO
      if (h.volumeId()%nro==0) continue;
      // reject all RO hit directly
      if (h.eDep()<0.0) continue;
      if (h.eDep()==0.0 && h.position()==CLHEP::Hep3Vector(0.0,0.0,0.0)) continue;
      hitmap[cal.getCrystalByRO(h.volumeId())].push_back(i);
    }

    CaloCrystalHitMCTruthCollection cr_hits;

    for(HitMap::const_iterator cr = hitmap.begin(); cr != hitmap.end(); ++cr) {

      // crystal id
      int cid = cr->first;

      // Prepare info for hit creation before the next iteration
      cr_hits.clear();

      // Loop over all hits found for this crystal
      vector<int> const& ihits = cr->second;

      for( size_t i=0; i<ihits.size(); i++ ) {
	StepPointMC const& h = (*hits)[ihits[i]];
	cr_hits.push_back(CaloCrystalHitMCTruth(cid,h.time(),h.eDep()));        
      }

      sort(cr_hits.begin(), cr_hits.end(), lessByTime<CaloCrystalHitMCTruthCollection::value_type>());

      // now form final hits if they are close enough in time

      CaloCrystalHitMCTruthCollection::value_type cHitMCTruth  = cr_hits[0];

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
DEFINE_FWK_MODULE(MakeCaloReadoutHits);
