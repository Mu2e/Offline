// See header file for description.

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/BFieldGeom/inc/BFMap.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/TrackerBFieldInfo.hh"

#include <iostream>

mu2e::TrackerBFieldInfo::TrackerBFieldInfo(){
  BFieldManager  const& bfm    = *GeomHandle<BFieldManager>();
  DetectorSystem const& detsys = *GeomHandle<DetectorSystem>();

  CLHEP::Hep3Vector trackerOriginWorld = detsys.toMu2e( CLHEP::Hep3Vector() );

  unsigned count{0};
  for ( auto m : bfm.getInnerMaps() ){
    if ( m->isValid( trackerOriginWorld ) ){
      _map = m.get();
      ++count;
    }
  }

  if ( count != 1 ){
    throw cet::exception("BField")
      << "TrackerBFieldInfo: could not locate exactly one BField innerMap containing the tracker origin. Number found: "
      << count
      << "\n Tracker origin in World: " << trackerOriginWorld;
  }

}
