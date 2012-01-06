#ifndef GeometryService_CosmicProductionPlane_hh
#define GeometryService_CosmicProductionPlane_hh
//
// Stores parameters for the production plane of the cosmic rays. These values are used by the cosmic ray generators. The WorldG4Maker uses theses values to check that no part of the production plane is outside of the world volume.
//
// $Id: CosmicProductionPlane.hh,v 1.1 2012/01/06 23:28:27 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/01/06 23:28:27 $
//
// Original authors Ralf Ehrlich and Rob Kutschke
//

// Mu2e includes.
#include "GeometryService/inc/Detector.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  // Forward references
  class SimpleConfig;

  class CosmicProductionPlane : public Detector
  {

  public:

    //sets default paramters which are used for the toy model
    CosmicProductionPlane();

    // Accept the compiler generator d'tor, copy c'tor and copy assigment.

    // implement Detector's method
    virtual std::string name() const { return "CosmicProductionPlane"; }

    double cosmicDx() const { return _cosmicDx; }
    double cosmicDz() const { return _cosmicDz; }
    double cosmicOffsetY() const { return _cosmicOffsetY; }

    //reads the config file (not the geo file) and sets the paramters for the DYB model
    void parametersDYB(SimpleConfig const& config) const;

  private:

    mutable double _cosmicDx, _cosmicDz, _cosmicOffsetY;
    //Since all detector components (like this object) can be accessed only via GeoHandle 
    //which returns only constant objects, all access functions need to be "const". 
    //However, these three variables can only be changed AFTER the detector has been constructed
    //(which will be done by cosmicDYB.cc). In order to allow them to be modified within a 
    //constant access function, the qualifier "mutable" was necessary.  
  };

} //namespace mu2e

#endif /* GeometryService_CosmicProductionPlane_hh */
