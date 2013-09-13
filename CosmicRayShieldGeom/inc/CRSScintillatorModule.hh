#ifndef CosmicRayShieldGeom_CRSScintillatorModule_hh
#define CosmicRayShieldGeom_CRSScintillatorModule_hh
//
// Representation of one Scintillator Module in CosmicRayShield
//

//
// $Id: CRSScintillatorModule.hh,v 1.7 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG somewhat based on Rob Kutschke' Sector
//

#include <vector>

#include "CosmicRayShieldGeom/inc/CRSScintillatorModuleId.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e 
{

  class CosmicRayShield;

  class CRSScintillatorModule
  {

    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSScintillatorModule();

    CRSScintillatorModule(CRSScintillatorModuleId const & id);

    // Accept the compiler generated destructor, copy constructor and
    // assignment operators

    const CRSScintillatorModuleId& id() const { return _id;}

    const std::vector<CRSScintillatorLayer>& getLayers() const
    {
      return _layers;
    }

    int nLayers() const
    {
      return _layers.size();
    }

    const CRSScintillatorLayer& getLayer ( int n ) const 
    {
      return _layers.at(n);
    }

    const CRSScintillatorLayer& getLayer ( CRSScintillatorLayerId const & lid) const 
    {
      return _layers.at(lid.getLayerNumber());
    }

    CRSScintillatorBar const & getBar ( const CRSScintillatorBarId& moduleid ) const
    {
      return _layers.at(moduleid.getLayerNumber()).getBar(moduleid);
    }

    // Formatted string embedding the id of the module.
    std::string name( std::string const & base ) const;

    // On readback from persistency, recursively recompute mutable members.
    //    void fillPointers ( const CosmicRayShield& cosmicRayShield ) const;

    private:

    CRSScintillatorModuleId _id;

    std::vector<CRSScintillatorLayer> _layers;
  };

}  //namespace mu2e
#endif /* CosmicRayShieldGeom_CRSScintillatorModule_hh */
