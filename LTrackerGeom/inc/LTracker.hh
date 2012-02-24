#ifndef LTrackerGeom_LTracker_hh
#define LTrackerGeom_LTracker_hh
//
// Hold all geometry and identifier information about
// an LTracker.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: LTracker.hh,v 1.13 2012/02/24 16:37:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:37:09 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "TrackerGeom/inc/Tracker.hh"

#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/SectorId.hh"

namespace mu2e {

  class LTracker: public Tracker{

    friend class LTrackerMaker;

  public:
    LTracker(){}
    ~LTracker(){}

    // Compiler generated copy and assignment constructors
    // should be OK.

    enum LTrackerDeviceId { undefined=-1, wedge, vane};

    double r0()   const { return _r0;}
    double z0()   const { return _z0;}
    double rOut() const { return _rOut;}

    double rInscribed() const { return _rInscribed;}

    double zHalfLength() const{
      return _halfLength;
    }

    double tiltY() const { return _tiltY;}
    double tiltX() const { return _tiltX;}

    // Depracated: kept for backwards compatibiltiy.
    double Tilt() const { return _tiltY;}

    int nSides() const { return _devices.at(0).getSectors().size(); }

    double strawRadius() const{
      return getStraw(StrawId(wedge,0,0,0)).getDetail().outerRadius();
    }

    std::string const& fillMaterial() const { return _fillMaterial; }

    // Check for legal identifiers.
    bool isLegal(DeviceId d) const{
      return ( d>-1 &&
               std::vector<Device>::size_type(d) <_devices.size()
               );
    }


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
    const std::vector<Device>& getDevices() const{
      return _devices;
    }

    const Device& getDevice ( DeviceId id) const{
      return _devices.at(id);
    }

    const Device& getSides() const {
      return _devices.at(wedge);
    }

    const Device& getVanes() const {
      return _devices.at(vane);
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
      return _strawDetail;
    }



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

    // Nominal values.
    // _r0 = Nominal radius of the center of the sector.
    // _z0 = position of the center of the tracker relative to the origin
    //       of the Mu2e coordinate system.
    // _Inscribed = inscribed radius that is tangent to the inner part of the sector.
    // tilts are described in Mu2e-doc-561-v2 ( or higher version).
    double _r0;
    double _z0;
    double _rInscribed;
    double _tiltY;
    double _tiltX;

    // Outer radius and half length ( in z ) of a logical volume that will
    // just contain the entire tracker.  Use to make the mother volume.
    double _rOut;
    double _halfLength;

    std::string _fillMaterial;

    // Detailed info about each type of straw.
    std::vector<StrawDetail> _strawDetail;

    // An LTracker is made of two devices, sides and vanes.
    std::vector<Device> _devices;

    // There will be pointers to the objects in this container.
    std::deque<Straw>  _allStraws;

    // Needed to complete the second phase of construction.
    void FillPointers1();
    void FillPointers2();

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers() const;

  };

} //namespace mu2e

#endif /* LTrackerGeom_LTracker_hh */
