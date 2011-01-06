#ifndef STRAW_HH
#define STRAW_HH
//
// Hold information about one straw in a tracker.
//
// $Id: Straw.hh,v 1.4 2011/01/06 22:37:11 wenzel Exp $
// $Author: wenzel $ 
// $Date: 2011/01/06 22:37:11 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "TrackerGeom/inc/StrawId.hh"
#include "TrackerGeom/inc/StrawIndex.hh"
#include "TrackerGeom/inc/StrawDetail.hh"
#include "TrackerGeom/inc/TubsParams.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

  // Forward declarations.
  class Tracker;


  class Straw{

    friend class Layer;
    friend class Sector;
    friend class Device;
    friend class LTracker;
    friend class TTracker;
    friend class LTrackerMaker;
    friend class TTrackerMaker;

  public:

    // A free function, returning void, that takes a const Straw& as an argument.
    typedef void (*StrawFunction)( const Straw& s);

    Straw();

    // Constructor using wire tangents.
    Straw( const StrawId& id,
           StrawIndex index,
           const CLHEP::Hep3Vector& c,
           const StrawDetail* detail,
           int    detailIndex,
           double wtx = 0.,
           double wty = 0.
           );

    // Constructor using wire unit vector.
    Straw( const StrawId& id,
           StrawIndex index,
           const CLHEP::Hep3Vector& c,
           const StrawDetail* detail,
           int detailIndex,
           CLHEP::Hep3Vector const& t
           );
  
    // Accept the compiler generated destructor, copy constructor and assignment operators

    const StrawId& Id() const { return _id;}
    StrawIndex Index() const { return _index;}

    const StrawDetail& getDetail() const { return *_detail;}

    // Return true if the argument is one of the nearest neighbours of this straw.
    bool isNearestNeighbour( StrawIndex idx ) const;

    const std::vector<StrawIndex>& nearestNeighboursByIndex() const{
      return _nearestByIndex;
    }
    const std::vector<StrawId>& nearestNeighboursById() const{
      return _nearestById;
    }
    // Formatted string embedding the id of the straw.
    std::string name( std::string const& base ) const;
  
    // Compiler generated copy and assignment constructors
    // should be OK.

    const CLHEP::Hep3Vector& getMidPoint() const {return _c;}

    // Delete one of these
    const CLHEP::Hep3Vector& getDirection() const { return _w;}
    const CLHEP::Hep3Vector& direction() const { return _w;}

    // Return G4TUBS parameters outer volume for this straw - gas volume.
    TubsParams getOuterTubsParams() const{
      return _detail->getOuterTubsParams();
    }

    // Return G4TUBS parameters for the straw skin.
    TubsParams getWallTubsParams() const{
      return _detail->getWallTubsParams();
    }
  
    // Return G4TUBS parameters for the wire.
    TubsParams getWireTubsParams() const{
      return _detail->getWireTubsParams();
    }

    // Straw Radius
    double getRadius() const {
      return _detail->outerRadius();
    }

    // Straw Thickness
    double getThickness() const {
      return _detail->thickness();
    }

    // Half length
    double getHalfLength() const {
      return _detail->halfLength();
    }

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

    int hack;

 protected:

    // Identifier
    StrawId _id;

    // Index into the array of all straws.
    StrawIndex _index;

    // Mid-point of the straw.
    CLHEP::Hep3Vector _c;

    // Detailed description of a straw.
    mutable const StrawDetail* _detail;
    int32_t _detailIndex;

    // Unit vector along the wire direction.
    // Need to add unit vectors along local u and v also.
    // Use Euler angle convention from G4.
    CLHEP::Hep3Vector _w;

    // Nearest neighbours.
    std::vector<StrawId>    _nearestById;
    std::vector<StrawIndex> _nearestByIndex;

  };

}  //namespace mu2e
#endif
