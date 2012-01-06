//
// Construct a CosmicProductionPlane object.
//
// $Id: CosmicProductionPlane.cc,v 1.1 2012/01/06 23:28:27 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/01/06 23:28:27 $
//
// Original authors Ralf Ehrlich
//

// Mu2e includes.
#include "GeometryService/inc/CosmicProductionPlane.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e
{
  CosmicProductionPlane::CosmicProductionPlane()
  {
    _cosmicDx=5000;
    _cosmicDz=5000;
    _cosmicOffsetY=0;
  }

  void CosmicProductionPlane::parametersDYB(SimpleConfig const& config) const
  {
    _cosmicDx=config.getDouble("cosmicDYB.dx",5000);
    _cosmicDz=config.getDouble("cosmicDYB.dz",5000);
    _cosmicOffsetY=config.getDouble("cosmicDYB.y0",0);
  }

}

