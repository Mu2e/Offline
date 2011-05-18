//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( LTracker or
// TTracker ). Return a pointer to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.cc,v 1.4 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e{

  const Tracker& getTrackerOrThrow(){

    art::ServiceHandle<GeometryService> geom;
    if( geom->hasElement<LTracker>() ){
      GeomHandle<LTracker> ltracker;
      return *ltracker;
    } else if ( geom->hasElement<TTracker>() ){
      GeomHandle<TTracker> ttracker;
      return *ttracker;
    } else if ( geom->hasElement<ITracker>() ){
      GeomHandle<ITracker> itracker;
      return *itracker;
    }


    throw cet::exception("GEOM")
      << "Expected one of L or T or I Trackers but found neither.\n";

  }

}
