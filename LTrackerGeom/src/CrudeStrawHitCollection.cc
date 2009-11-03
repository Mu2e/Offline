#include <iostream>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"

// Mu2e includes
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/CrudeStrawHitCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

using namespace std;

namespace mu2e{

  CrudeStrawHitCollection::
  CrudeStrawHitCollection( edm::Event const& event,
			   edm::Handle<CrudeStrawHitPData> const& hits ):
    _event(&event),
    _hits(&(*hits)),
    _index()
  {
    FillIndex();
  }

  CrudeStrawHitCollection::
  CrudeStrawHitCollection( edm::Event const& event,
			   CrudeStrawHitPData const& hits ):
    _event(&event),
    _hits(&hits),
    _index()
  {
    FillIndex();
  }

  // Fill the _index member datum.
  void CrudeStrawHitCollection::FillIndex(){

    // Geometry for the LTracker.
    GeomHandle<LTracker> ltracker;
    int nChannels = ltracker->getAllStraws().size();
    
    // Fill the fast indexing array.
    _index.assign( nChannels, -1 );
    for ( int i=0; i<_hits->size(); ++i){
      CrudeStrawHit const& hit = (*_hits)[i];
      StrawIndex idx = hit.strawIndex;
      _index[idx.asInt()] = i;
    }
  }

  void CrudeStrawHitCollection::getStepPointMC( int i,
						vector<StepPointMC const*>& v ) const{

    // The requested hit.
    CrudeStrawHit const& hit = _hits->at(i);

    // In future, change this to check to see if an indirect path
    // to the StepPointMC is available.  If so, navigate through
    // that path.
    if ( hit.precursorType != CrudeStrawHit::stepPointMC){
      throw cms::Exception("Hits")
	<< "Requested precursor for a CrudeStrawHit and that precursor could not be found.\n"
	<< "Requested type: StepPointMC "
	<< "Available type code: "
	<< hit.precursorType;
    }

    // Fill the return argument.
    hit.getStepPointMC( *_event, v);

  }

  int CrudeStrawHitCollection::strawIndexToHitIndexOrThrow( StrawIndex idx ) const{

    int ihit(_index[idx.asInt()]);

    if ( ihit == -1 ) {
      throw cms::Exception("RANGE")
	<< "There is no hit for this StrawIndex: "
	<< idx
	<< "\n";
    }

    return ihit;
  }
  

}

