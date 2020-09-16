/*------------------------------------------------------------

A filter module that do not write events with no
hits in tracker and calorimeter.
A integer parameter defined in the py file can exclude
one of the detectors from the filter.
0 skip all events with no hit in both detector
1 skip only events with no hits in the tracker
2 skip events with no hit in the calorimeter


Original author Giovanni Onorato


-----------------------------------------------------------*/

// C++ includes
#include <iostream>
#include <stdexcept>
#include <string>

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
//#include <boost/shared_ptr.hpp>
#include "fhiclcpp/ParameterSet.h"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
        
  class FilterEmptyEvents : public art::EDFilter {
  public:
      explicit FilterEmptyEvents(fhicl::ParameterSet const& pset);
//      virtual ~FilterEmptyEvents();
      virtual bool filter(art::Event& e ) override;

  private:
    // Control parameter: 0 to filter both tracker and calorimeter
    //                    1 to filter only tracker
    //                    2 to filter only calorimeter

    int _diagLevel;
    int _keepTrackOrCalo;
    string _generatorModuleLabel;
    string _makerModuleLabel;
    string _caloReadoutModuleLabel;
    unsigned _minTHits, _minCHits;

    bool _hasTHits;
    bool _hasCHits;
  };
                  
  //constructor
  FilterEmptyEvents::FilterEmptyEvents(fhicl::ParameterSet const& pset)
    : EDFilter{pset},
    _diagLevel(pset.get<int>("diagLevel",0)),
    _keepTrackOrCalo(pset.get<int>("keepTrackOrCalo",1)),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
    _minTHits(pset.get<unsigned>("MinTrackerHits",0)),
    _minCHits(pset.get<unsigned>("MinCaloHits",0))
    {}
                  
  bool FilterEmptyEvents::filter(art::Event& e ) {

    if (_keepTrackOrCalo < 0 || _keepTrackOrCalo>2) {
      cout << "Meaningless KeepTrackOrCalo parameter value."
           << " It must be 0, 1 or 2. We use 0 for this case "
           << "as default value." << endl;
    }

    //Get Tracker hits and set a boolean true if there are hits.
    art::Handle<StrawHitCollection> pdataHandle;
    e.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    _hasTHits = (hits->size()>_minTHits);

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    e.getByLabel(_generatorModuleLabel, genParticles);

    if (!_hasTHits) {
      if (_diagLevel > 0) {
        cout << "Filtered Event " << e.id().event() << " has " << hits->size() << " hits in the tracker," << endl;
      }
    }

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if(geom->hasElement<DiskCalorimeter>()) {

      // Get handles to calorimeter collections
      art::Handle<CaloHitCollection> caloHits;
      e.getByLabel(_caloReadoutModuleLabel, caloHits);

      //Set a boolean true if there are calorimeter hits
      _hasCHits = (caloHits.isValid() && (caloHits->size() > _minCHits));

      if (!_hasCHits) {
        if (_diagLevel > 0) {
          cout << "Filtered Event " << e.id().event() << " has " << caloHits->size() << " hits in the calorimeter" << endl;
        }
      }
    } else {
      //Set a boolean false if there is no Calorimeter
      _hasCHits = false;
//      cout << "No calorimeter in the geometry" << endl;
    }

    if (_keepTrackOrCalo == 1) {
      return (_hasTHits);
    }

    if (_keepTrackOrCalo == 2) {
      return (_hasCHits);
    }

    if (_keepTrackOrCalo == 0) {
      return (_hasTHits&&_hasCHits);
    }

    return false;

  }

}

DEFINE_ART_MODULE(mu2e::FilterEmptyEvents);
