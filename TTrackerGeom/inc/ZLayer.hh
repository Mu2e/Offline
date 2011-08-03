#ifndef TrackerGeom_ZLayer_hh
#define TrackerGeom_ZLayer_hh
//
// Holds information about one layer (ordered by Z) in a tracker.
//

//
// $Id: ZLayer.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <vector>
#include <iostream>

#include "TTrackerGeom/inc/ZLayerId.hh"
#include "TrackerGeom/inc/Straw.hh"


namespace mu2e {

  class TTracker;
  
  class ZLayer{

    friend class TTracker;
    friend class TTrackerMaker;

  public:

    ZLayer():_id(ZLayerId()) {}
    ZLayer( const ZLayerId& id, double z ): _id(id), _z(z) {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const ZLayerId& id() const { return _id;}

    // The next layer down in the hierarchy is Straws
    int   nStraws()                           const { return _straws.size(); }
    const std::vector<Straw const *>& getStraws () const { return _straws; }
    Straw const* getStrawptr ( int n)         const { return _straws.at(n); }
    Straw const* getStrawptr  ( const StrawId& stid ) const{
      return _straws.at(stid.getStraw());
    }

    // Formatted string embedding the id of the ZLayer.
    std::string name( std::string const& base ) const;

  // Get geometric abtraction information about this ZLayer
 
  double midZ() const {return _z;};

   // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& Ttracker ) const;

  protected:

    ZLayerId _id;
    double _z;
    std::vector<Straw const *> _straws;

 };

  inline
  std::ostream & operator<< (std::ostream & os, ZLayer const &x) 
    { return os << x.name("ZLayer "); }

}  //namespace mu2e
#endif /* TrackerGeom_ZLayer_hh */
