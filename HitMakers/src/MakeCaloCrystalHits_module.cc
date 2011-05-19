//
// An EDProducer Module that reads CaloHit objects and turns them into
// CaloCrystalHit objects, collection
//
// $Id: MakeCaloCrystalHits_module.cc,v 1.4 2011/05/19 23:51:50 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/19 23:51:50 $
//
// Original author KLG
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

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
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"

#include "Mu2eUtilities/inc/sort_functors.hh"

using namespace std;
using art::Event;

namespace mu2e {

  class MakeCaloCrystalHits : public art::EDProducer {

  public:

    explicit MakeCaloCrystalHits(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _maximumEnergy(pset.get<double>("maximumEnergy",1000.0)), //MeV
      _minimumTimeGap(pset.get<double>("minimumTimeGap",100.0)),// ns
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _messageCategory("CaloHitMaker")

    {
      // Tell the framework what we make.
      produces<CaloCrystalHitCollection>();
    }
    virtual ~MakeCaloCrystalHits() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Parameters
    double _minimumEnergy;  // minimum energy in the RO to count it

    double _maximumEnergy;  // energy of a saturated RO

    double _minimumTimeGap; // to merge the hits

    string _g4ModuleLabel;  // Name of the module that made the input hits.

    string _caloReadoutModuleLabel; // Name of the module that made the calo hits.

    // A category for the error logger.
    const std::string _messageCategory;

    void fixEnergy(CaloCrystalHitCollection::value_type & caloCrystalHit,
                   int tnro, double electronEdep);

  };

  void MakeCaloCrystalHits::beginJob(){
  }

  void MakeCaloCrystalHits::fixEnergy(CaloCrystalHitCollection::value_type & caloCrystalHit,
                                      int tnro, double electronEdep) {

    int nridu = caloCrystalHit.numberOfROIdsUsed();

    if ( _diagLevel > 0 ) {
      cout << __func__ << ": fixing energy: " << caloCrystalHit.energyDep()
           << ", used roids: " << nridu
           << ", energyDepT: " << caloCrystalHit.energyDepTotal() << endl;
    }

    if (nridu > 0 && nridu < tnro ) {
      caloCrystalHit.setEnergyDep(caloCrystalHit.energyDep()/
                                  double(nridu)*double(tnro));
    }

    // fix only if all ro are saturated
    if (nridu == 0 && caloCrystalHit.energyDepTotal()/tnro >= electronEdep) {
      caloCrystalHit.setEnergyDep(caloCrystalHit.energyDepTotal());
    }



    if ( _diagLevel > 0 ) {
      cout << __func__ << ": fixed  energy: " <<  caloCrystalHit.energyDep()
           << ", used roids: " << nridu
           << ", energyDepT: " << caloCrystalHit.energyDepTotal() << endl;
    }

    return;

  }

  void MakeCaloCrystalHits::produce(art::Event& event) {

    // Source of the info
    art::Handle<CaloHitCollection> caloHits;
    // A container to hold the output hits.
    auto_ptr<CaloCrystalHitCollection> caloCrystalHits(new CaloCrystalHitCollection);

    if ( _diagLevel > 0 ) cout << __func__ << ": begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Get handles to calorimeter (RO aka APD) collection

    event.getByLabel(_caloReadoutModuleLabel, caloHits);
    bool haveCalo = caloHits.isValid();

    if ( _diagLevel > -1 && !haveCalo) cout << __func__ << ": No CaloHits" << endl;

    if( !haveCalo) return;

    if (caloHits->size()<=0) {
      if ( _diagLevel > 0 ) cout << __func__ << ": 0 CaloHits" << endl;
      // Add the empty hit collection to the event
      event.put(caloCrystalHits);
      return;
    }

    // Get calorimeter geometry description
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) {
      throw cet::exception("GEOM")
        << "Expected calorimeter, but found none";
    }
    //    GeomHandle<Calorimeter> cg;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    int nro = cal.nROPerCrystal();
    double electronEdep = cal.getElectronEdep();

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      _diagLevel > 0 &&
        cout << __func__ << ": Total number of hit RO = " << caloHits->size() << endl;
      for( size_t i=0; i<caloHits->size(); ++i ) {
        _diagLevel > 0 &&
          cout << __func__ << ": " << (*caloHits)[i]
               << " CrystalId: " << cal.getCrystalByRO((*caloHits)[i].id()) << endl;
      }
    }

    // hits should be organized by crystal
    // should they be also separated by time?

    // Product Id of the input points.
    art::ProductID const& caloHitCollId(caloHits.id());
    //mcptr.push_back(DPIndex(id,straw_hits[0]._hit_id));

    // Instatiate caloHitsSorted from caloHits which is const,
    //  we may do it with pointers to hits later instead
    CaloHitCollection caloHitsSorted(*caloHits);

    // Sort them by id & time
    sort(caloHitsSorted.begin(), caloHitsSorted.end(),
         lessByCIdAndTime<CaloHitCollection::value_type>(cal));

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << __func__ << ": Total number of hit RO Sorted = " << caloHitsSorted.size() << endl;
      for(CaloHitCollection::const_iterator i = caloHitsSorted.begin();
          i != caloHitsSorted.end(); ++i) {
        //      for( size_t i=0; i<caloHitsSorted.size(); ++i ) {
        cout << __func__ << ": ";
        i->print(cout,false);
        cout << " CrystalId: " << cal.getCrystalByRO(i->id()) << endl;
      }
    }

    // generate the CaloCrystalHits
    // collect same crystal id/time hits

    // we will use all the hits, but the energyDep will be collected only for the good ones
    // energyDepTotal will include all the energy

    // the energy has to be between min and maximum

    CaloHitCollection::value_type const & hit0 = *caloHitsSorted.begin();
    _diagLevel > 0 && cout << __func__ << ": Original RO hit:  " << hit0 <<endl;

    CaloCrystalHitCollection::value_type caloCrystalHit;

    if ( hit0.energyDep()>= _minimumEnergy && hit0.energyDep() < _maximumEnergy ) {
      caloCrystalHit.assign(cal.getCrystalByRO(hit0.id()),caloHitCollId,hit0);
    } else {
      caloCrystalHit.assignEnergyToTot(cal.getCrystalByRO(hit0.id()),caloHitCollId,hit0);
    }

    _diagLevel > 0 && cout << __func__ << ": As in CaloCrystalHit: "
                           << caloCrystalHit << endl;

    for( CaloHitCollection::const_iterator i = caloHitsSorted.begin()+1;
         i != caloHitsSorted.end(); ++i) {

      CaloHitCollection::value_type const & hit = *i;

      _diagLevel > 0 && cout << __func__ << ": Original RO hit: " << hit << endl;

      int cid = cal.getCrystalByRO(hit.id());

      _diagLevel > 0 &&
        cout << __func__ << ": old, new cid:  " << caloCrystalHit.id() << ", " << cid << endl;
      _diagLevel > 0 &&
        cout << __func__ << ": old, new time: " << caloCrystalHit.time() << ", " << hit.time() << endl;

      _diagLevel > 0 &&
        cout << __func__ << ": time difference, gap: "
             << (hit.time() - caloCrystalHit.time()) << ", "
             << _minimumTimeGap << endl;

      if (caloCrystalHit.id() == cid &&
          (( hit.time() - caloCrystalHit.time()) < _minimumTimeGap) ) {

        // here we decide if the hit is "good"

        if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) {

          caloCrystalHit.add(caloHitCollId, hit);

        } else {

          caloCrystalHit.addEnergyToTot(caloHitCollId, hit);

        }

        _diagLevel > 0 && cout << __func__ << ": Added to the hit:  " << caloCrystalHit << endl;

      } else {


        fixEnergy(caloCrystalHit,nro,electronEdep);

        if (caloCrystalHit.energyDep()>0.0) {
          _diagLevel > 0 && cout << __func__ << ": Inserting old hit: " << caloCrystalHit << endl;
          (*caloCrystalHits).push_back(caloCrystalHit);
        }

        // this resets the caloCrystalHit and sets its id and puts one hit in

        if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) {

          caloCrystalHit.assign(cid, caloHitCollId, hit);

        } else {

          caloCrystalHit.assignEnergyToTot(cid, caloHitCollId, hit);

        }
        _diagLevel > 0 && cout << __func__ << ": Created new hit:   " << caloCrystalHit << endl;

      }

    }

    if (_diagLevel > 1) {
      cout << __func__ << ": roIds of last old hit:";
      cout << " size: " << caloCrystalHit.roIds().size();
      for( size_t i=0; i<caloCrystalHit.roIds().size(); ++i) {
        cout  << " " << i << ": " << caloCrystalHit.roIds()[i];
      }
      cout << endl;
    }

    fixEnergy(caloCrystalHit,nro,electronEdep);

    if (caloCrystalHit.energyDep()>0.0) {
      _diagLevel > 0 && cout << __func__ << ": Inserting last old hit: " << caloCrystalHit << endl;
      (*caloCrystalHits).push_back(caloCrystalHit);
    }

    if ( _diagLevel > 0 ) cout << __func__ << ": (*caloCrystalHits).size() "
                               << (*caloCrystalHits).size() << endl;

    // Add the output hit collection to the event (invalidates caloCrystalHits handle)
    event.put(caloCrystalHits);

    if ( _diagLevel > 0 ) cout << __func__ << ": ncalls " << ncalls << endl;
    if ( _diagLevel > 0 ) cout << __func__ << ": end" << endl;
    return;

  } // end of ::produce

}

using mu2e::MakeCaloCrystalHits;
DEFINE_ART_MODULE(MakeCaloCrystalHits);
