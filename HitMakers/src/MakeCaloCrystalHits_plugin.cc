//
// An EDProducer Module that reads CaloHit objects and turns them into
// CaloCrystalHit objects, collection
//
// $Id: MakeCaloCrystalHits_plugin.cc,v 1.1 2010/10/28 20:43:58 genser Exp $
// $Author: genser $
// $Date: 2010/10/28 20:43:58 $
//
// Original author KLG
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
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"

#include "ToyDP/inc/CaloCrystalHitCollection.hh"
//#include "ToyDP/inc/CaloCrystalHitMCTruthCollection.hh"


using namespace std;
using edm::Event;

namespace mu2e {

  // utility functor to sort hits by time

  class lessByIdAndTimeCaloHits {

  public:
    
    bool operator() (const CaloHit& a, const CaloHit& b) const {
      return (a.roId() < b.roId() ||
              (a.roId() == b.roId() &&
               a.time() < b.time() 
               ) 
              );
    }

  };

  class MakeCaloCrystalHits : public edm::EDProducer {
  public:
    explicit MakeCaloCrystalHits(edm::ParameterSet const& pset) : 

      // Parameters

      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _minimumEnergy(pset.getUntrackedParameter<double>("minimumEnergy",0.0001)), // MeV
      _minimumTimeGap(pset.getUntrackedParameter<double>("minimumTimeGap",100.0)),// ns
      _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
      _messageCategory("CaloHitMaker")

    {
      // Tell the framework what we make.
      produces<CaloCrystalHitCollection>();
      //      produces<CaloCrystalHitMCTruthCollection>();
    }
    virtual ~MakeCaloCrystalHits() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void produce( edm::Event& e, edm::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Parameters
    double _minimumEnergy;  // minimum energy in the APD

    double _minimumTimeGap;
      
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void MakeCaloCrystalHits::beginJob(edm::EventSetup const& ){
  }

  void MakeCaloCrystalHits::produce(edm::Event& event, edm::EventSetup const&) {

    if ( _diagLevel > 0 ) cout << "MakeCaloCrystalHits: produce() begin" << endl;
      
    static int ncalls(0);
    ++ncalls;

    // Get handles to calorimeter (RO aka APD) collections
    edm::Handle<CaloHitCollection>        caloHits;
    //    edm::Handle<CaloHitMCTruthCollection> caloHitsMCTruth;
    event.getByType(caloHits);
    //    event.getByType(caloMC);
    //    bool haveCalo = ( caloHits.isValid() && caloMC.isValid() );
    bool haveCalo = caloHits.isValid();

    if ( _diagLevel > -1 && !haveCalo) cout << "MakeCaloCrystalHits: No CaloHits" << endl;

    if( !haveCalo) return;

    // A container to hold the output hits.
    auto_ptr<CaloCrystalHitCollection>        caloCrystalHits(new CaloCrystalHitCollection);

    if (caloHits->size()<=0) {
      if ( _diagLevel > 0 ) cout << "MakeCaloCrystalHits: 0 CaloHits" << endl;
      // Add the empty hit collection to the event
      event.put(caloCrystalHits);
      return;
    }

    // Get calorimeter geometry description
    edm::Service<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) {
      throw cms::Exception("GEOM")
        << "Expected calorimeter, but found none";
    }
    GeomHandle<Calorimeter> cg;
    //    int nro = cg->nROPerCrystal();

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      _diagLevel > 0 && 
        cout << "MakeCaloCrystalHits: Total number of hit RO = " << caloHits->size() << endl;
      for( size_t i=0; i<caloHits->size(); ++i ) {
        _diagLevel > 0 && 
          cout << "MakeCaloCrystalHits: " << (*caloHits)[i]
               << " CrystalId: " << cg->getCrystalByRO((*caloHits)[i].roId()) << endl;
      }
    }

    // hits should be organized by crystal
    // should they be also separated by time?

    //    auto_ptr<CaloCrystalHitMCTruthCollection> truthHits(new CaloCrystalHitMCTruthCollection);

    // Product Id of the input points.
    edm::ProductID const& caloHitCollId(caloHits.id());
    //mcptr.push_back(DPIndex(id,straw_hits[0]._hit_id));

    // Instatiate caloHitsTO from caloHits which is const, we may do it with pointers to hits later instead
    CaloHitCollection caloHitsTO(*caloHits);

    // Sort them by id & time
    sort(caloHitsTO.begin(), caloHitsTO.end(), lessByIdAndTimeCaloHits() );

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeCaloCrystalHits: Total number of hit RO TO = " << caloHitsTO.size() << endl;
      for(CaloHitCollection::const_iterator i = caloHitsTO.begin(); i != caloHitsTO.end(); ++i) {
        //      for( size_t i=0; i<caloHitsTO.size(); ++i ) {
        cout << "MakeCaloCrystalHits: ";
        i->print(cout,false);
        cout << " CrystalId: " << cg->getCrystalByRO(i->roId()) << endl;
      }
    }

    // generate the CaloCrystalHits
    // collect same time/crystal id hits

    CaloHit const & hit0 = caloHitsTO[0];
    _diagLevel > 0 && cout << "MakeCaloCrystalHits: Original APD hit:  " << hit0 <<endl;

    CaloCrystalHit caloCrystalHit(cg->getCrystalByRO(hit0.roId()), caloHitCollId, hit0);
    
    _diagLevel > 0 && cout << "MakeCaloCrystalHits: As in CaloCrystalHit: "
                           << caloCrystalHit << endl;

    for( size_t i=1; i<caloHitsTO.size(); ++i) {
      
      CaloHit const & hit = caloHitsTO.at(i);

      _diagLevel > 0 && cout << "MakeCaloCrystalHits: Original APD hit: " << hit << endl;
      
      int cid = cg->getCrystalByRO(hit.roId());

      _diagLevel > 0 && 
        cout << "MakeCaloCrystalHits: old, new cid:  " << caloCrystalHit.crystalId() << ", " << cid << endl;
      _diagLevel > 0 && 
        cout << "MakeCaloCrystalHits: old, new time: " << caloCrystalHit.time() << ", " << hit.time() << endl;

      _diagLevel > 0 && 
        cout << "MakeCaloCrystalHits: time difference, gap: " 
             << (hit.time() - caloCrystalHit.time()) << ", "
             << _minimumTimeGap << endl;

      if (caloCrystalHit.crystalId() == cid && 
          (( hit.time() - caloCrystalHit.time()) < _minimumTimeGap) ) {

        caloCrystalHit.add(caloHitCollId, hit);
        _diagLevel > 0 && cout << "MakeCaloCrystalHits: Added to the hit:  " << caloCrystalHit << endl;

      } else {

        _diagLevel > 0 && cout << "MakeCaloCrystalHits: Inserting old hit: " << caloCrystalHit << endl;

        (*caloCrystalHits).push_back(caloCrystalHit);
        // this resets the caloCrystalHit and sets its id and puts one hit in
        caloCrystalHit.assign(cid, caloHitCollId, hit);
        _diagLevel > 0 && cout << "MakeCaloCrystalHits: Created new hit:   " << caloCrystalHit << endl;

      }

    }

    if (_diagLevel > 1) {
      cout << "MakeCaloCrystalHits: roIds of last old hit:";
      cout << " size: " << caloCrystalHit.roIds().size();
      for( size_t i=0; i<caloCrystalHit.roIds().size(); ++i) {
        cout  << " " << i << ": " << caloCrystalHit.roIds()[i];
      }
      cout << endl;
    }

    _diagLevel > 0 && cout << "MakeCaloCrystalHits: Inserting last old hit: " << caloCrystalHit << endl;

    (*caloCrystalHits).push_back(caloCrystalHit);
  
    if ( _diagLevel > 0 ) cout << "MakeCaloCrystalHits: (*caloCrystalHits).size() " 
                               << (*caloCrystalHits).size() << endl;

    // Add the output hit collection to the event (invalidates caloCrystalHits handle)
    event.put(caloCrystalHits);

    if ( _diagLevel > 0 ) cout << "MakeCaloCrystalHits: ncalls " << ncalls << endl;
    if ( _diagLevel > 0 ) cout << "MakeCaloCrystalHits: produce() end" << endl;
    return;

  } // end of ::analyze.
  
}

using mu2e::MakeCaloCrystalHits;
DEFINE_FWK_MODULE(MakeCaloCrystalHits);
