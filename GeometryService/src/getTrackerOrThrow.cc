//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( TTracker or
// ITracker ). Return a pointer to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.cc,v 1.5 2013/03/26 23:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/26 23:28:23 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e{

  const Tracker& getTrackerOrThrow(){

    art::ServiceHandle<GeometryService> geom;
    if ( geom->hasElement<TTracker>() ){
      GeomHandle<TTracker> ttracker;
      return *ttracker;
    } else if ( geom->hasElement<ITracker>() ){
      GeomHandle<ITracker> itracker;
      return *itracker;
    }

    throw cet::exception("GEOM")
      << "Expected one of T or I Tracker but found neither.\n";

  }

}
