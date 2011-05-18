/*------------------------------------------------------------

A filter module that do not write events with no
hits in tracker and calorimeter.
A integer parameter defined in the py file can exclude
one of the detectors from the filter.
0 skip all events with no hit in both detector
1 skip only events with no hits in the tracker
2 skip events with no hit in the calorimeter

$Id: FilterEmptyEvents_module.cc,v 1.3 2011/05/18 02:27:14 wb Exp $
$Author: wb $
$Date: 2011/05/18 02:27:14 $

Original author Giovanni Onorato


-----------------------------------------------------------*/

// C++ includes
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <string>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
//#include <boost/shared_ptr.hpp>
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"


using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class FilterEmptyEvents : public art::EDFilter {
  public:
    explicit FilterEmptyEvents(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _keepTrackOrCalo(pset.get<int>("keepTrackOrCalo",1)),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")){
    }
    virtual ~FilterEmptyEvents() {
    }
    virtual bool filter(art::Event& e );

  private:

    // Control parameter: 0 to filter both tracker and calorimeter
    //                    1 to filter only tracker
    //                    2 to filter only calorimeter

    int _diagLevel;
    int _keepTrackOrCalo;
    string _makerModuleLabel;

    bool _hasTHits;
    bool _hasCHits;

  };

  bool FilterEmptyEvents::filter(art::Event& e ) {

    if (_keepTrackOrCalo < 0 || _keepTrackOrCalo>2) {
      cout << "Meaningless KeepTrackOrCalo parameter value."
           << " It must be 0, 1 or 2. We use 0 for this case "
           << "as default value." << endl;
    }

    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    // const Tracker& tracker = getTrackerOrThrow();

    //Get Tracker hits and set a boolean true if there are hits.
    art::Handle<StrawHitCollection> pdataHandle;
    e.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    _hasTHits = (hits->size()>0);

    // Get handles to the generated and simulated particles.
    art::Handle<ToyGenParticleCollection> genParticles;
    e.getByType(genParticles);

    if (!_hasTHits) {
      if (_diagLevel > 0) {
        cout << "Event " << e.id().event() << " has no hits in the tracker" << endl;
      }
    }

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( geom->hasElement<Calorimeter>() ) {

      // Get handles to calorimeter collections
      art::Handle<CaloHitCollection> caloHits;
      e.getByType(caloHits);

      //Set a boolean true if there are calorimeter hits
      _hasCHits = (caloHits.isValid() && (caloHits->size() > 0));

      if (!_hasCHits) {
        if (_diagLevel > 0) {
          cout << "Event " << e.id().event() << " has no hits in the calorimeter" << endl;
        }
      }
    } else {
      //Set a boolean false if there is no Calorimeter
      _hasCHits = false;
      cout << "No calorimeter in the geometry" << endl;
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

using mu2e::FilterEmptyEvents;
DEFINE_ART_MODULE(FilterEmptyEvents);
