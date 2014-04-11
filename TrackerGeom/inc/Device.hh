#ifndef TrackerGeom_Device_hh
#define TrackerGeom_Device_hh
//
// Hold information about one device in a tracker.
//
// $Id: Device.hh,v 1.10 2014/04/11 04:39:13 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:39:13 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>

// Mu2e includes
#include "TrackerGeom/inc/DeviceId.hh"
#include "TrackerGeom/inc/Sector.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Tracker;

  class Device{

    friend class TTracker;
    friend class TTrackerMaker;

  public:

    // A free function, returning void, that takes a const Device& as an argument.
    typedef void (*DeviceFunction)( const Device& s);

    Device():_id(-1),_rotation(0.),_origin(),_sectors(),_exists(true){}

    explicit Device( const DeviceId& id,
            CLHEP::Hep3Vector const& origin = CLHEP::Hep3Vector(0.,0.,0.),
            double rotation = 0., bool exists = true ):
      _id(id),
      _rotation(rotation),
      _origin(origin),
      _exists(exists) {
    }

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    DeviceId id() const { return _id;}

    double rotation() const { return _rotation; }

    const CLHEP::Hep3Vector & origin() const { return _origin; }

    int nSectors() const{
      return _sectors.size();
    }

    const std::vector<Sector>& getSectors () const{
      return _sectors;
    }

    const Sector& getSector ( int n) const {
      return _sectors.at(n);
    }

    const Sector& getSector ( const SectorId& sid ) const{
      return _sectors.at(sid.getSector());
    }

    const Layer& getLayer ( const LayerId& lid ) const{
      return _sectors.at(lid.getSector()).getLayer(lid);
    }

    const Straw& getStraw ( const StrawId& sid ) const{
      return _sectors.at(sid.getSector()).getStraw(sid);
    }

    // Formatted string embedding the id of the sector.
    std::string name( std::string const& base ) const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

    bool exists() const {
      return _exists;
    }

#ifndef __CINT__

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
            i !=e; ++i){
        i->forAllStraws(f);
      }
    }

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllLayers ( F& f) const{
      for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
            i !=e; ++i){
        i->forAllLayers(f);
      }
    }

    template <class F>
    inline void forAllSectors ( F& f) const{
      for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
            i !=e; ++i){
        f(*i);
      }
    }

#endif

  protected:

    DeviceId            _id;
    double              _rotation;
    CLHEP::Hep3Vector   _origin;
    std::vector<Sector> _sectors;
    bool                _exists;
  };

} //namespace mu2e

#endif /* TrackerGeom_Device_hh */
