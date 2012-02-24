#ifndef TTrackerGeom_TTracker_hh
#define TTrackerGeom_TTracker_hh
//
// Hold all geometry and identifier information about
// a TTracker.  This is intended as a "data only"
// class.
//
// $Id: TTracker.hh,v 1.11 2012/02/24 16:37:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:37:09 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "TTrackerGeom/inc/Manifold.hh"
#include "TTrackerGeom/inc/Support.hh"

#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/SectorId.hh"
#include "TrackerGeom/inc/StrawDetail.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/TubsParams.hh"

#include "TTrackerGeom/inc/Station.hh"

namespace mu2e {

  class TTracker: public Tracker{

    friend class TTrackerMaker;

  public:
    TTracker(){}  // TODO: insert proper initializer list, starting w/ base class

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    void fillPointers () const;

    double rOut() const { return _rOut;}
    double z0()   const { return _z0;}
    double zHalfLength() const;

    double strawRadius() const{
      return getStraw(StrawId(0,0,0,0)).getDetail().outerRadius();
    }

    std::string const& envelopeMaterial() const { return _envelopeMaterial; }

    // Check for legal identifiers.
    bool isLegal(DeviceId d) const{
      return ( d>-1 &&
               std::vector<Device>::size_type(d) <_devices.size()
               );
    };

    bool isLegal(const SectorId& sid) const{
      return (isLegal(sid.getDeviceId()) &&
              sid.getSector() >-1   &&
              std::vector<Sector>::size_type(sid.getSector()) < getDevice(sid.getDeviceId()).getSectors().size()
              );
    }

    typedef std::vector<Sector>::size_type stypeLayer;
    bool isLegal(const LayerId& lid ) const{
      return ( isLegal(lid.getSectorId()) &&
               lid.getLayer() > -1   &&
               std::vector<Layer>::size_type(lid.getLayer()) < getSector(lid.getSectorId()).getLayers().size()
               );
    }

    bool isLegal(const StrawId& sid) const{
      return ( isLegal(sid.getLayerId()) &&
               sid.getStraw() > -1       &&
               std::vector<Straw>::size_type(sid.getStraw()) < getLayer(sid.getLayerId()).getStraws().size()
               );
    }

    // Accessors
    int nDevices() const{
      return _devices.size();
    }

    const std::vector<Device>& getDevices() const{
      return _devices;
    }

    const Device& getDevice ( DeviceId id) const{
      return _devices.at(id);
    }

    const Sector& getSector ( const SectorId& sid ) const{
      return _devices.at(sid.getDevice()).getSector(sid);
    }

    const Layer& getLayer ( const LayerId& lid ) const{
      return _devices.at(lid.getDevice()).getLayer(lid);
    }

    const Straw& getStraw ( const StrawId& sid ) const{
      return _devices.at(sid.getDevice()).getStraw(sid);
    }

    const Straw& getStraw ( StrawIndex i ) const{
      return _allStraws.at(i.asInt());
    }

    int nStations() const{
      return _stations.size();
    }

    const std::vector<Station>& getStations() const{
      return _stations;
    }

    const Station& getStation ( StationId id) const{
      return _stations.at(id);
    }

    const std::deque<Straw>& getAllStraws() const {return _allStraws;}

    const std::vector<StrawDetail>& getStrawDetails() const{
      return _strawDetails;
    }

    const Support& getSupportParams () const{
      return _supportParams;
    }

    const std::vector<double>& getManifoldHalfLengths () const{
      return _manifoldHalfLengths;
    }

    // Shape parameters of an envelop holding one device.
    TubsParams getDeviceEnvelopeParams() const;

    TubsParams getTrackerEnvelopeParams() const;


#ifndef __CINT__

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
            i !=e; ++i){
        i->forAllStraws(f);
      }
    }

    template <class F>
    inline void forAllLayers ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
            i !=e; ++i){
        i->forAllLayers(f);
      }
    }

    template <class F>
    inline void forAllSectors ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
            i !=e; ++i){
        i->forAllSectors(f);
      }
    }

    template <class F>
    inline void forAllDevices ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
            i !=e; ++i){
        f(*i);
      }
    }

#endif


  protected:

    // Position of the center of the tracker, in the Mu2e coordinate system.
    double _z0;

    // Outer radius of a logical volume that will just contain the entire tracker.
    double _rOut;

    std::string _envelopeMaterial;

    // Detailed info about each type of straw.
    std::vector<StrawDetail> _strawDetails;

    // An TTracker is made of two devices, sides and vanes.
    std::vector<Device> _devices;

    // An alternative viewpoint:  
    // A TTracker is made of a collection of Stations.
    std::vector<Station> _stations;

    // There will be pointers to the objects in this container.
    std::deque<Straw>  _allStraws;

    std::vector<Manifold> _allManifolds;

    // All supports are the same shape.
    Support _supportParams;

    // All manifolds are the same shape.
    std::vector<double> _manifoldHalfLengths;

    // Inner radius of inside edge of innermost straw.
    double _envelopeInnerRadius;

  };

} //namespace mu2e

#endif /* TTrackerGeom_TTracker_hh */
