//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( TTracker or
// ITracker ). Return a pointer to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.cc,v 1.6 2014/09/03 16:39:15 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/03 16:39:15 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"

namespace mu2e{

  const Tracker& getTrackerOrThrow(){

    art::ServiceHandle<GeometryService> geom;
    if ( geom->hasElement<TTracker>() ){
      GeomHandle<TTracker> ttracker;
      return *ttracker;
    }

    throw cet::exception("GEOM")
      << "Expected T Tracker but didn't find neither.\n";

  }

}
