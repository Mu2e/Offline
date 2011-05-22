#ifndef CosmicRayShieldGeom_CRSScintillatorShield_hh
#define CosmicRayShieldGeom_CRSScintillatorShield_hh
//
// Representation of one ScintillatorShield in CosmicRayShield.
//
// $Id: CRSScintillatorShield.hh,v 1.4 2011/05/22 20:28:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/22 20:28:13 $
//
// Original author KLG based on Rob Kutschke's Device
//

// C++ includes
#include <vector>

// Mu2e includes
#include "CosmicRayShieldGeom/inc/CRSScintillatorShieldId.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CRSScintillatorShield{

    friend class CosmicRayShieldMaker;

  public:

    CRSScintillatorShield():_id(-1){}

    CRSScintillatorShield(
                          CRSScintillatorShieldId const & id,
                          std::string       const & name,
                          CLHEP::Hep3Vector const & localOffset,  // offset in Hall Air
                          std::vector<double> const & globalRotationAngles,
                          CLHEP::Hep3Vector const & globalOffset, // offset in World
                          double const              halfThickness,
                          std::vector<int>  const & numberOfModules);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const CRSScintillatorShieldId id() const { return _id;}

    //    const double rotation() const { return _rotation; }

    CLHEP::Hep3Vector const & getLocalOffset() const { return _localOffset; }

    std::vector<double> const & getGlobalRotationAngles() const { return _globalRotationAngles;}

    CLHEP::Hep3Vector const & getGlobalOffset() const { return _globalOffset; }

    int nModules() const{
      return _numberOfFullModules+_numberOfHalfModules;
    }

    int nCRSScintillatorFullModules() const{
      return _numberOfFullModules;
    }

    const std::vector<CRSScintillatorModule>& getCRSScintillatorModules() const{
      return _modules;
    }

    const CRSScintillatorModule& getModule( int n) const {
      return _modules.at(n);
    }

    const CRSScintillatorModule& getModule( const CRSScintillatorModuleId& moduleid ) const{
      return _modules.at(moduleid.getModuleNumber());
    }

    const CRSScintillatorLayer& getLayer( const CRSScintillatorLayerId& lid ) const{
      return _modules.at(lid.getModuleNumber()).getLayer(lid);
    }

    CRSScintillatorBar const & getBar( const CRSScintillatorBarId& bid ) const{
      return _modules.at(bid.getModuleNumber()).getBar(bid);
    }

    // Formatted string embedding the id of the shield.

    std::string name( std::string const & base ) const;

    // On readback from persistency, recursively recompute mutable members.
    //    void fillPointers( const CosmicRayShield& cosmicRayShield ) const;

  private:

    CRSScintillatorShieldId _id;

    std::string _name;

    // position in the parent frame
    CLHEP::Hep3Vector _localOffset;

    std::vector<double> _globalRotationAngles;

    // position in Mu2e
    CLHEP::Hep3Vector _globalOffset;

    // outer dimensions; the thickness for now
    double _halfThickness;

    // there are Full and Half Modules, each module "knows it"
    // we need to know how many there will be ahead of adding them

    int _numberOfFullModules;
    int _numberOfHalfModules;

    std::vector<CRSScintillatorModule> _modules;

  };

} //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorShield_hh */
