#ifndef MU2EHALL_HH
#define MU2EHALL_HH

//====================================================================================
//
// Geometry of the hall, dirt, etc.
//
// Kyle Knoepfel, 2014-15
//
// This class stores all geometries relevant for constructing the Mu2e
// hall and dirt volumes.  The geometries are constructed in several
// steps using the Geometry service:
//
//   (1) Mu2eHallMaker::makeBuilding(...) constructs all concrete volumes that
//       comprise the building.
//
//   (2) The three-dimensional envelope of the Mu2e building is then
//       determined from the minimum/maximum extents of the building volumes.
//
//       ***NB***: the envelope is that of the CONCRETE volumes, and
//                 does not include any dirt volumes.
//
//   (3) The x, y, and z Mu2e envelope values are used to determine
//       the boundaries for the dirt volumes, which are constructed
//       using Mu2eHallMaker::makeDirt(...).
//
//   (4) The dirt grade level in the y direction is calculated
//       whenever the world is created, using the
//       "dirtDsAreaFirstFloorS" solid.  This is determined in the
//       body of WorldG4Maker::make(...), which is guaranteed to be
//       called after the dirt volumes have already been loaded.
//
//   (5) The building and dirt volumes are constructed in G4 using
//       Mu2eG4/src/constructHall.cc
//
//====================================================================================

#include <map>
//#include <vector>

//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/TwoVector.h"

#include "Offline/GeomPrimitives/inc/ExtrudedSolid.hh"
#include "Offline/GeomPrimitives/inc/RotExtrudedSolid.hh"
#include "Offline/GeomPrimitives/inc/GenericTrap.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class Mu2eHallMaker;

  class Mu2eHall : virtual public Detector {
  public:

    const std::map<std::string,ExtrudedSolid>& getBldgSolids()      const { return bldgSolids_; }
    const std::map<std::string,ExtrudedSolid>& getDirtSolids()      const { return dirtSolids_; }
    const std::map<std::string,RotExtrudedSolid>& getRotSolids()      const { return rotatedSolids_; }
    const std::map<std::string,GenericTrap>&   getDirtTrapSolids()  const { return dirtTrapSolids_; }

    const ExtrudedSolid&
    getBldgSolid( const std::string& str ) const;

    const ExtrudedSolid&
    getDirtSolid( const std::string& str ) const;

    const RotExtrudedSolid&
    getRotSolid( const std::string& str ) const;

    const GenericTrap&
    getDirtTrapSolid( const std::string& str ) const;

    double getWallExtentz( const std::string& , const int  ) const;

    //----------------------------------------------------------------
  private:
    friend class Mu2eHallMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eHall(){}

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    std::map<std::string,ExtrudedSolid> bldgSolids_;
    std::map<std::string,ExtrudedSolid> dirtSolids_;
    std::map<std::string,RotExtrudedSolid> rotatedSolids_;
    std::map<std::string,GenericTrap>   dirtTrapSolids_;

  };

}

#endif/*MU2EHALL_HH*/
