#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Hold all geometry and identifier information about
// a Tracker.  This is intended as a "data only" class.
//  Note that the only 'aligned' information is implicit in the individual straws
//
// An un-aligned version is provided by GeometryService
// and an aligned verison is provided  by ProditionsService
//
// Original author Rob Kutschke
//

#include <vector>
#include <array>
#include <limits>
#include <string>
#include <memory>

#include "cetlib_except/exception.h"

#include "Mu2eInterfaces/inc/Detector.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

#include "TrackerGeom/inc/Plane.hh"
#include "TrackerGeom/inc/Panel.hh"
#include "TrackerGeom/inc/StrawProperties.hh"
#include "TrackerGeom/inc/G4Tracker.hh"

namespace mu2e {
  class Tracker : public Detector, public ProditionsEntity {
    using xyzVec = CLHEP::Hep3Vector; // switch to root XYZVec TODO
    public:
    typedef std::shared_ptr<Tracker> ptr_t;
    typedef std::shared_ptr<const Tracker> cptr_t;
    using PlaneCollection = std::array<Plane,StrawId::_nplanes>;
    using PanelCollection = std::array<Panel,StrawId::_nupanels>;
    using StrawCollection = std::array<Straw,StrawId::_nustraws>;
    using StrawIndexMap = std::array<uint16_t,StrawId::_maxval>;

    using PEType = std::array<bool,StrawId::_nplanes>;

    // default constructor results in non-functional object, but is required by proditions service

    constexpr static const char* cxname = {"Tracker"};
    // copy constructor
    Tracker(const Tracker& other);
    // construct from a set of straws and their global properties.  It would be better to take ownership of
    // existing straws via std::move, but the need for the copy constructor above precludes that
    Tracker(StrawCollection const& straws, StrawProperties const& sprops,const std::shared_ptr<G4Tracker>& g4tracker,
    PEType const& pexists);

    // accessors
    // origin in nominal tracker coordinate system
    const xyzVec& origin() const { return _origin; }

    // common straw attributes that don't depend on the specific element
    const StrawProperties& strawProperties() const { return _strawprops; }

    // constituent access
    size_t nPlanes() const { return _planes.size(); }
    size_t nPanels() const { return _panels.size(); }
    size_t nStraws() const { return _straws.size(); }
    PlaneCollection const& planes() const { return _planes; }
    PanelCollection const& panels() const{ return _panels; }
    StrawCollection const& straws() const{ return _straws; }
    // fast indexed lookup by StrawId through indirect map
    uint16_t strawIndex(const StrawId& id) const { return _strawindex.at(id.asUint16()); }
    const Plane& plane( const StrawId& id ) const{ return _planes.at(id.getPlane()); }
    const Panel& panel( const StrawId& id ) const{ return _panels.at(id.uniquePanel()); }
    const Straw& straw( const StrawId& id) const{ return _straws[strawIndex(id)]; }

    // access the G4Tracker
    const std::shared_ptr<G4Tracker>& g4Tracker() const { return _g4tracker; }

    // deprecated interface: do not write new code using these, and replace existing calls opportunistically
    const Plane& getPlane( const StrawId& id ) const{ return _planes.at(id.getPlane()); }
    const Plane& getPlane( uint16_t n ) const{ return _planes.at(n); }
    PlaneCollection const& getPlanes() const { return _planes; }
    const Panel& getPanel( const StrawId& id ) const{ return _panels.at(id.uniquePanel()); }
    const Straw& getStraw( const StrawId& id) const{ return _straws[strawIndex(id)]; }
    Straw& getStraw( const StrawId& id) { return _straws[strawIndex(id)]; }
    StrawCollection const& getStraws() const{ return _straws; }
    // the following are deprecated: access should be through StrawProperties
    double strawInnerRadius() const{ return _strawprops._strawInnerRadius; }
    double strawOuterRadius() const{ return _strawprops._strawOuterRadius; }
    double strawWallThickness() const{ return _strawprops._strawWallThickness; }
    double outerMetalThickness() const{ return _strawprops._outerMetalThickness; }
    double innerMetal1Thickness() const{ return _strawprops._innerMetal1Thickness; }
    double innerMetal2Thickness() const{ return _strawprops._innerMetal2Thickness; }
    double wireRadius()           const { return _strawprops._wireRadius; }
    double wirePlateThickness()   const { return _strawprops._wirePlateThickness; }

    // depracted 'exists' interface: should switch to TrackerStatus FIXME!
    auto const& planesExist() const { return _planeExists; }
    bool planeExists(StrawId const& id) const { return _planeExists[id.plane()]; }

    private:
    xyzVec _origin;
    // global straw properties
    StrawProperties _strawprops;
    // Dense arrays
    PlaneCollection _planes;
    PanelCollection _panels;
    // fundamental geometric content is in the following
    StrawCollection _straws;
    // indirection from StrawId, for efficient lookup
    StrawIndexMap _strawindex;
    // plane existence: use cases of this should switch to using TrackerStatus and this should be removed FIXME!!
    PEType _planeExists;
    // g4 content
    std::shared_ptr<G4Tracker> _g4tracker;
  };

} //namespace mu2e

#endif /* TrackerGeom_Tracker_hh */
