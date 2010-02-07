#ifndef DEVICE_HH
#define DEVICE_HH

//
// Hold information about one device in a tracker.
//

//
// $Id: Device.hh,v 1.1 2010/02/07 00:29:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/07 00:29:41 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "TrackerGeom/inc/DeviceId.hh"
#include "TrackerGeom/inc/Sector.hh"

namespace mu2e {

class Device{

  friend class LTracker;
  friend class LTrackerMaker;

public:

  // A free function, returning void, that takes a const Device& as an argument.
  typedef void (*DeviceFunction)( const Device& s);

  Device():_id(-1){}
  Device( const DeviceId& id ):_id(id){}
  ~Device(){}
 
  // Compiler generated copy and assignment constructors
  // should be OK.
  
  const DeviceId Id() const { return _id;}

  const std::vector<Sector>& getSectors () const{ 
    return _sectors;
  }

  const Sector& getSector ( int n) const { 
    return _sectors.at(n);
  }

  const Sector& getSector ( const SectorId& sid ) const{
    return _sectors.at(sid._sector);
  }

  const Layer& getLayer ( const LayerId& lid ) const{
    return _sectors.at(lid.getSector()).getLayer(lid);
  }

  const Straw& getStraw ( const StrawId& sid ) const{
    return _sectors.at(sid.getSector()).getStraw(sid);
  }

  // Formatted string embedding the id of the sector.
  std::string name( std::string const& base ) const;


#ifndef __CINT__

  // Loop over all straws and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void Device::forAllStraws ( F& f) const{
    for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
	  i !=e; ++i){
      i->forAllStraws(f);
    }
  }

  // Loop over all straws and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void Device::forAllLayers ( F& f) const{
    for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
	  i !=e; ++i){
      i->forAllLayers(f);
    }
  }

  template <class F>
  inline void Device::forAllSectors ( F& f) const{
    for ( std::vector<Sector>::const_iterator i=_sectors.begin(), e=_sectors.end();
	  i !=e; ++i){
      f(*i);
    }
  }

#endif

protected:

  DeviceId _id;
  std::vector<Sector> _sectors;

};

} //namespace mu2e

#endif
