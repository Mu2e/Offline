#ifndef TrackerGeom_Straw_hh
#define TrackerGeom_Straw_hh
//
// Hold information about one straw in a tracker.
//
// Original author Rob Kutschke
//

#include <vector>

#include "DataProducts/inc/StrawId.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Tracker;

  class Straw{

    friend class TrackerMaker;
    friend class AlignedTrackerMaker;

  public:

    // A free function, returning void, that takes a const Straw& as an argument.
    typedef void (*StrawFunction)( const Straw& s);

    Straw();

    // Constructor using wire tangents.
    Straw( const StrawId& id,
           const CLHEP::Hep3Vector& c,
           double wtx = 0.,
           double wty = 0.
           );

    // Constructor using wire unit vector.
    Straw( const StrawId& id,
           const CLHEP::Hep3Vector& c,
           CLHEP::Hep3Vector const& t
           );

    // Accept the compiler copy constructor and assignment operators

    // I don't think that this class should have virtual functions but
    // since it does, it must have a virtual destructor.
    virtual ~Straw(){}


    const StrawId& id() const { return _id;}

    // Navigation, normally within a panel:

    // The following routnies will return the nearest neighbors (or the
    // single nearest neighbor if this straw is at the end of its layer):
    const std::vector<StrawId>& nearestNeighboursById() const{
      return _nearestById;
    }
    const std::vector<StrawId>& preampNeighboursById() const{
      return _preampById;
    }

    // Formatted string embedding the id of the straw.
    std::string name( std::string const& base ) const;

    virtual const CLHEP::Hep3Vector& getMidPoint() const {return _c;}

    // Delete one of these
    virtual  const CLHEP::Hep3Vector& getDirection() const { return _w;}
    const CLHEP::Hep3Vector& direction() const { return _w;}

    // outer radius
    virtual double getRadius() const;

    virtual double innerRadius() const;

    virtual double  getThickness() const;

    // Half length
    virtual double halfLength() const;

    // active length is a little smaller
    virtual double activeHalfLength() const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker* tracker ) const;

    int hack;
    bool operator==(const Straw other) const {
      return _id == other.id();
    }
    bool operator>(const Straw other) const {
      return _id > other.id();
    }
   bool operator<(const Straw other) const {
      return _id < other.id();
   }
 protected:

    // Identifier
    StrawId _id;

    // Mid-point of the straw.
    CLHEP::Hep3Vector _c;

    // Unit vector along the wire direction.
    // Need to add unit vectors along local u and v also.
    // Use Euler angle convention from G4.
    CLHEP::Hep3Vector _w;

    mutable const Tracker* _tracker;

    // Nearest neighbours.
    std::vector<StrawId>    _nearestById;
    std::vector<StrawId>    _preampById;

  };

}  //namespace mu2e
#endif /* TrackerGeom_Straw_hh */
