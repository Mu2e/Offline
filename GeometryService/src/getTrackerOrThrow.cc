//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( LTracker or
// TTracker ). Return a pointer to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.cc,v 1.1 2010/04/18 00:32:42 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/04/18 00:32:42 $
//
// Original author Rob Kutschke
//

// Framework includes 
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"

namespace mu2e{

  const Tracker& getTrackerOrThrow(){

    edm::Service<GeometryService> geom;
    if( geom->hasElement<LTracker>() ){
      GeomHandle<LTracker> ltracker;
      return *ltracker;
    } else if ( geom->hasElement<TTracker>() ){
      GeomHandle<TTracker> ttracker;
      return *ttracker;
    }

    throw cms::Exception("GEOM") 
      << "Expected one of L or T Trackers but found neither.\n";

  }

}
