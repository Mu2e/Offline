#ifndef STRAW_HH
#define STRAW_HH
//
// Hold information about one straw in a tracker.
//
//
// $Id: Straw.hh,v 1.3 2009/10/28 13:40:46 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/28 13:40:46 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "LTrackerGeom/inc/StrawId.hh"
#include "LTrackerGeom/inc/StrawIndex.hh"
#include "LTrackerGeom/inc/StrawDetail.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

class Straw{

  friend class Layer;
  friend class Sector;
  friend class Device;
  friend class LTracker;
  friend class LTrackerMaker;


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
  
  ~Straw ();

  const StrawId& Id() const { return _id;}
  StrawIndex Index() const { return _index;}

  const StrawDetail& getDetail() const { return *_detail;}

  const std::vector<const Straw*>& nearestNeighbours() const{
    return _nearest;
  }

  // Return true if the argument is one of the nearest neighbours of this straw.
  bool isNearestNeighbour( StrawIndex idx ) const;

  const std::vector<StrawIndex>& nearestNeighboursByIndex() const{
    return _nearestByIndex;
  }
  
  // Compiler generated copy and assignment constructors
  // should be OK.

  const CLHEP::Hep3Vector& getMidPoint() const {return _c;}

  const CLHEP::Hep3Vector& getDirection() const { return _w;}

  int hack;
  
protected:

  // Identifier
  StrawId _id;

  // Index into the array of all straws.
  StrawIndex _index;

  // Mid-point of the straw.
  CLHEP::Hep3Vector _c;

  // Detailed description of a straw.
  const StrawDetail* _detail;
  int _detailIndex;

  // Unit vector along the wire direction.
  // Need to add unit vectors along local u and v also.
  // Use Euler angle convention from G4.
  CLHEP::Hep3Vector _w;

  // Nearest neighbours.
  std::vector<const Straw *> _nearest;

  std::vector<StrawId> _nearestById;
  std::vector<StrawIndex> _nearestByIndex;

  // Vector of all straws, from the enclosing Tracker class.
  //static std::deque<Straw>* _allStraws;

};

}  //namespace mu2e

#endif
