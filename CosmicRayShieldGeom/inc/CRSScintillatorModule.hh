#ifndef CosmicRayShieldGeom_CRSScintillatorModule_hh
#define CosmicRayShieldGeom_CRSScintillatorModule_hh
//
// Representation of one Scintillator Module in CosmicRayShield
//

//
//
// Original author KLG somewhat based on Rob Kutschke' Sector
//

#include <vector>

#include "CosmicRayShieldGeom/inc/CRSScintillatorModuleId.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CRSAbsorberLayer.hh"
#include "CosmicRayShieldGeom/inc/CRSAluminumSheet.hh"
#include "CosmicRayShieldGeom/inc/CRSFEB.hh"


namespace mu2e 
{

  class CosmicRayShield;

  class CRSScintillatorModule
  {

    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    CRSScintillatorModule();

    public:

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

    const std::vector<CRSAbsorberLayer>& getAbsorberLayers() const
    {
      return _absorberLayers;
    }

    int nAbsorberLayers() const
    {
      return _absorberLayers.size();
    }

    const CRSAbsorberLayer& getAbsorberLayer ( int n ) const 
    {
      return _absorberLayers.at(n);
    }

    const std::vector<CRSAluminumSheet>& getAluminumSheets() const
    {
      return _aluminumSheets;
    }

    int nAluminumSheets() const
    {
      return _aluminumSheets.size();
    }

    const CRSAluminumSheet& getAluminumSheet ( int n ) const 
    {
      return _aluminumSheets.at(n);
    }

    const std::vector<CRSFEB>& getFEBs() const
    {
      return _FEBs;
    }

    int nFEBs() const
    {
      return _FEBs.size();
    }

    const CRSFEB& getFEB ( int n ) const 
    {
      return _FEBs.at(n);
    }


    // Formatted string embedding the id of the module.
    std::string name( std::string const & base ) const;

    // On readback from persistency, recursively recompute mutable members.
    //    void fillPointers ( const CosmicRayShield& cosmicRayShield ) const;

    private:

    CRSScintillatorModuleId _id;

    std::vector<CRSScintillatorLayer> _layers;
    std::vector<CRSAbsorberLayer> _absorberLayers;
    std::vector<CRSAluminumSheet> _aluminumSheets;
    std::vector<CRSFEB> _FEBs;
  };

}  //namespace mu2e
#endif /* CosmicRayShieldGeom_CRSScintillatorModule_hh */
