#ifndef TrackerGeom_Straw_hh
#define TrackerGeom_Straw_hh
//
// Hold information about one straw in a tracker.
//
// $Id: Straw.hh,v 1.18 2013/03/26 23:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/26 23:28:23 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "TrackerGeom/inc/StrawDetail.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward declarations.
  class Tracker;


  class Straw{

    friend class Layer;
    friend class Panel;
    friend class Plane;
    friend class TTracker;
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

    // Accept the compiler copy constructor and assignment operators

    // I don't think that this class should have virtual functions but 
    // since it does, it must have a virtual destructor.
    virtual ~Straw(){}


    const StrawId& id() const { return _id;}
    StrawIndex index() const { return _index;}
    int detailIndex() const { return _detailIndex;}

    const StrawDetail& getDetail() const { return *_detail;}

    // Navigation, normally within a panel:

    // Return true if the argument is one of the nearest neighbours of this straw.
    bool isNearestNeighbour( StrawIndex idx ) const;

    // The following routnies will return the nearest neighbors (or the
    // single nearest neighbor if this straw is at the end of its layer):
    const std::vector<StrawIndex>& nearestNeighboursByIndex() const{
      return _nearestByIndex;
    }
    const std::vector<StrawId>& nearestNeighboursById() const{
      return _nearestById;
    }
    
    
    // The following routines deal with straws in this straw's layer
    // as well as the other layer in the same panel.  They return 
    // the enum NO_STRAW when the appropriate requested straw does not exist.
    StrawIndex nextOuterSameLayer() const {
      return  _nextOuterL;
    }
    StrawIndex nextInnerSameLayer() const {
      return  _nextInnerL;
    }
    StrawIndex nextOuterInPanel() const {
      return _nextOuterP;  // This will always be in the opposite layer
    }
    StrawIndex nextInnerInPanel() const {
      return _nextInnerP;
    }

    // Formatted string embedding the id of the straw.
    std::string name( std::string const& base ) const;

    virtual const CLHEP::Hep3Vector& getMidPoint() const {return _c;}

    // Delete one of these
    virtual  const CLHEP::Hep3Vector& getDirection() const { return _w;}
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
    virtual double getRadius() const {
      return _detail->outerRadius();
    }

    // Straw Thickness
    virtual double getThickness() const {
      return _detail->thickness();
    }

    // Half length
    virtual double getHalfLength() const {
      return _detail->halfLength();
    }

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

    int hack;
    bool operator==(const Straw other) const {
      return _index == other.index();
    }
    bool operator>(const Straw other) const {
      return _index > other.index();
    }
   bool operator<(const Straw other) const {
      return _index < other.index();
   }
 protected:

    // Identifier
    StrawId _id;

    // Index into the array of all straws.
    StrawIndex _index;

    // Mid-point of the straw.
    CLHEP::Hep3Vector _c;

    // Detailed description of a straw.
    mutable const StrawDetail* _detail;
    int _detailIndex;

    // Unit vector along the wire direction.
    // Need to add unit vectors along local u and v also.
    // Use Euler angle convention from G4.
    CLHEP::Hep3Vector _w;

    // Nearest neighbours.
    std::vector<StrawId>    _nearestById;
    std::vector<StrawIndex> _nearestByIndex;

    StrawIndex _nextOuterL;
    StrawIndex _nextInnerL;
    StrawIndex _nextOuterP;
    StrawIndex _nextInnerP;
  };

}  //namespace mu2e
#endif /* TrackerGeom_Straw_hh */
