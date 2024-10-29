#ifndef CosmicRayShieldGeom_CosmicRayShield_hh
#define CosmicRayShieldGeom_CosmicRayShield_hh

//
// Representation of CosmicRayShield
//
//
// Original author KLG
//

// c++ includes
#include <string>

// clhep includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Includes from Mu2e
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSSupportStructure.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"


namespace mu2e
{

  // Forward reference.
  class SimpleConfig;
  class CosmicRayShieldMaker;

  class CosmicRayShield : virtual public Detector
  {

    friend class CosmicRayShieldMaker;

    public:

    CosmicRayShield() {}

    // Get ScintillatorShield
    CRSScintillatorShield const & getCRSScintillatorShield(const CRSScintillatorShieldId& id) const
    {
      return _scintillatorShields.at(id);
    }

    CRSScintillatorModule const & getModule( const CRSScintillatorModuleId& moduleid ) const
    {
      return _scintillatorShields.at(moduleid.getShieldNumber()).getModule(moduleid);
    }

    CRSScintillatorLayer const & getLayer( const CRSScintillatorLayerId& lid ) const
    {
      return _scintillatorShields.at(lid.getShieldNumber()).getLayer(lid);
    }

    CRSScintillatorBar const & getBar( const CRSScintillatorBarId& bid ) const
    {
      return _scintillatorShields.at(bid.getShieldNumber()).getBar(bid);
    }

    std::vector<CRSScintillatorShield> const & getCRSScintillatorShields() const
    {
      return _scintillatorShields;
    }

    std::vector<std::shared_ptr<CRSScintillatorBar> > const & getAllCRSScintillatorBars() const
    {
      return _allCRSScintillatorBars;
    }

    const CRSScintillatorBar& getBar ( CRSScintillatorBarIndex index ) const
    {
      return *_allCRSScintillatorBars.at(index.asInt());
    }

    //sector names are e.g. R1 (for only the R1 sector), or R (for all R sectors)
    std::vector<double> getSectorHalfLengths(const std::string &sectorName) const;
    CLHEP::Hep3Vector getSectorPosition(const std::string &sectorName) const;

    std::vector<CRSSupportStructure> const & getSupportStructures() const {return _supportStructures;}

    const std::string& getName() const { return _name; }

    private:

    std::string _name;  // name of loaded geometry to coordinate with conditions

    std::vector<CRSScintillatorShield>                _scintillatorShields;    //Every "shield" holds a vector of modules.
                                                                               //Every module hold a vector of layers.
                                                                               //Every layer holds a vector of pointers
                                                                               //to CRV bars.
    std::vector<std::shared_ptr<CRSScintillatorBar> > _allCRSScintillatorBars; //This vector holds pointers to all CRV bars,
                                                                               //(the same objects used in all layers).

    void getMinMaxPoints(const std::string &sectorName, std::vector<double> &minPoint, std::vector<double> &maxPoint) const;

    std::vector<CRSSupportStructure> _supportStructures;

  };

}
#endif /* CosmicRayShieldGeom_CosmicRayShield_hh */
