#ifndef GeometryService_src_CosmicProductionPlaneMaker_hh
#define GeometryService_src_CosmicProductionPlaneMaker_hh
//
// Construct a CosmicProductionPlane object.
//
// $Id: CosmicProductionPlaneMaker.hh,v 1.1 2012/01/06 23:28:27 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/01/06 23:28:27 $
//
// Original authors Ralf Ehrlich and Rob Kutschke
//

// C++ includes
#include <memory>

// Mu2e includes.
#include "GeometryService/inc/Detector.hh"
#include "GeometryService/inc/CosmicProductionPlane.hh"

namespace mu2e {

  class CosmicProductionPlaneMaker: public Detector
  {

  public:

    CosmicProductionPlaneMaker() 
    {
      _cosmicProductionPlane = std::auto_ptr<CosmicProductionPlane>(new CosmicProductionPlane());
    }

    // Accept the compiler generator d'tor, copy c'tor and copy assigment.

    // Give control of the cosmic production plane to the caller.
    std::auto_ptr<CosmicProductionPlane> getCosmicProductionPlanePtr() 
    { 
      return _cosmicProductionPlane;
    }

  private:

    std::auto_ptr<CosmicProductionPlane> _cosmicProductionPlane;

  };

} //namespace mu2e

#endif /* GeometryService_src_CosmicProductionPlaneMaker_hh */
