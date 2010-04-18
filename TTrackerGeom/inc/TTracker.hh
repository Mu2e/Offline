#ifndef TTRACKER_HH
#define TTRACKER_HH
//
// Hold all geometry and identifier information about
// a TTracker.  This is intended as a "data only"
// class.
//
// $Id: TTracker.hh,v 1.1 2010/04/18 00:37:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:37:16 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <deque>

#include "TTrackerGeom/inc/Support.hh"
#include "TTrackerGeom/inc/Manifold.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/StrawDetail.hh"
#include "TrackerGeom/inc/SectorId.hh"
#include "TrackerGeom/inc/TubsParams.hh"


namespace mu2e {

  class TTracker: public Tracker{

  friend class TTrackerMaker;

  public:
    TTracker(){}
    ~TTracker(){};
    
    // Compiler generated copy and assignment constructors
    // should be OK.

    void fillPointers () const;

    virtual std::string name() const { return "TTracker";}

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
      return (isLegal(sid._did) && 
	      sid._sector >-1   &&
	      std::vector<Sector>::size_type(sid._sector) < getDevice(sid._did).getSectors().size()
	      );
    }

    typedef std::vector<Sector>::size_type stypeLayer;
    bool isLegal(const LayerId& lid ) const{
      return ( isLegal(lid._sid) &&
	       lid._layer > -1   &&
	       std::vector<Layer>::size_type(lid._layer) < getSector(lid._sid).getLayers().size()
	       );
    }

    bool isLegal(const StrawId& sid) const{
      return ( isLegal(sid._lid) &&
	       sid._n > -1       &&
	       std::vector<Straw>::size_type(sid._n) < getLayer(sid._lid).getStraws().size()
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

    const std::deque<Straw>& getAllStraws() const {return _allStraws;}

    const std::vector<StrawDetail>& getStrawDetails() const{
      return _strawDetails;
    }

    // Shape parameters of an envelop holding one device.
    TubsParams getDeviceEnvelopeParams() const;

    TubsParams getTrackerEnvelopeParams() const;


#ifndef __CINT__

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void TTracker::forAllStraws ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	    i !=e; ++i){
	i->forAllStraws(f);
      }
    }

    template <class F>
    inline void TTracker::forAllLayers ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	    i !=e; ++i){
	i->forAllLayers(f);
      }
    }

    template <class F>
    inline void TTracker::forAllSectors ( F& f) const{
      for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	    i !=e; ++i){
	i->forAllSectors(f);
      }
    }
    
    template <class F>
    inline void TTracker::forAllDevices ( F& f) const{
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

#endif
