#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Hold all geometry and identifier information about
// a Tracker.  This is intended as a "data only" class.
//
// An un-aligned version can be made by GeometryService
// and an aligned verison can be made by ProditionsService
// by copying then adding alignment to the GeometryService version
//
// Original author Rob Kutschke
//

#include <vector>
#include <array>
#include <limits>
#include <string>

#include "cetlib_except/exception.h"


#include "Mu2eInterfaces/inc/Detector.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "TrackerGeom/inc/Manifold.hh"
#include "TrackerGeom/inc/Support.hh"
#include "TrackerGeom/inc/SupportModel.hh"
#include "TrackerGeom/inc/SupportStructure.hh"

#include "TrackerGeom/inc/Plane.hh"
#include "DataProducts/inc/PanelId.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PlacedTubs.hh"

namespace mu2e {

  class Tracker : public Detector, public ProditionsEntity {

        friend class TrackerMaker;

  public:

    typedef std::shared_ptr<Tracker> ptr_t;
    typedef std::shared_ptr<const Tracker> cptr_t;

    constexpr static uint16_t _maxRedirect =
      ((StrawId::_nplanes -1) << StrawId::_planesft) +
      ((StrawId::_npanels -1) << StrawId::_panelsft) +
      StrawId::_nstraws;

    // =============== NewTracker Public Objects End   ==============

    Tracker():_name("AlignedTracker") {
      if (StrawId::_nlayers != 2)
        throw cet::exception("GEOM")
          << "Expect configuration with 2 layers per panel\n";
    }  // TODO: insert proper initializer list, starting w/ base class

    // Need to be able to copy, but this has internal
    // pointers, so needs a deep copy, but can except default dtor

    // copy constructor
    Tracker(const Tracker& other);

    std::string const& name() const { return _name; }

    void fillPointers () const;

    double rOut() const { return _rOut;}
    double z0()   const { return _z0;}
    double zHalfLength() const;

    double strawInnerRadius() const{ return _strawInnerRadius; }
    double strawOuterRadius() const{ return _strawOuterRadius; }
    double strawWallThickness() const{ return _strawWallThickness; }
    double outerMetalThickness() const{ return _outerMetalThickness; }
    double innerMetal1Thickness() const{ return _innerMetal1Thickness; }
    double innerMetal2Thickness() const{ return _innerMetal2Thickness; }
    double wireRadius()           const { return _wireRadius; }
    double wirePlateThickness()   const { return _wirePlateThickness; }

    std::string const& wallMaterialName()            const{ return _wallMaterialName; }
    std::string const& wallCoreMaterialName()        const{ return  wallMaterialName();  }
    std::string const& wallOuterMetalMaterialName()  const{ return _outerMetalMaterial;  }
    std::string const& wallInnerMetal1MaterialName() const{ return _innerMetal1Material; }
    std::string const& wallInnerMetal2MaterialName() const{ return _innerMetal2Material; }
    std::string const& gasMaterialName()             const{ return _gasMaterialName; }
    std::string const& wireMaterialName()            const{ return _wireMaterialName; }
    std::string const& wireCoreMaterialName()        const{ return  wireMaterialName();  }
    std::string const& wirePlateMaterialName()       const{ return _wirePlateMaterial;   }

    // istraw is StrawId::straw()
    double getStrawHalfLength(int istraw) const { return _strawHalfLengths[istraw];}
    double getStrawActiveHalfLength(int istraw) const { return _strawActiveHalfLengths[istraw];}

    std::string const& envelopeMaterial() const { return _envelopeMaterial; }

    // Check for legal identifiers. (for what decays to StrawId)
    bool isLegal(const StrawId strid) const{
       return strid.valid();
    }

    // Accessors
    uint16_t nPlanes() const { return StrawId::_nplanes; }

    uint16_t nPanels() const { return StrawId::_nupanels; }

    uint16_t nStraws() const { return StrawId::_nustraws; }

    const std::vector<double>& getManifoldHalfLengths () const{
      return _manifoldHalfLengths;
    }

    double panelOffset() const { return _panelZOffset; }

    SupportModel getSupportModel() const{
      return _supportModel;
    }

    const Support& getSupportParams () const{
      return _supportParams;
    }

    const SupportStructure& getSupportStructure() const{
      return _supportStructure;
    }


    TubsParams getPlaneEnvelopeParams() const{
      return _planeEnvelopeParams;
    }

    TubsParams getPanelEnvelopeParams() const{
      return _panelEnvelopeParams;
    }

    const TubsParams& getInnerTrackerEnvelopeParams() const{
      return _innerTrackerEnvelopeParams;
    }

    PlacedTubs mother() const{
      return _mother;
    }


    const Plane& getPlane( const StrawId id ) const{
      return _planes.at(id.getPlane());
    }
    const Plane& getPlane( uint16_t n ) const{
      return _planes.at(n);
    }

    std::array<Plane,StrawId::_nplanes> const& getPlanes() const {
      return _planes;
    }

    const Panel& getPanel( const StrawId id ) const{
      return _panels.at(id.uniquePanel());
    }

    const Straw& getStraw( const StrawId id) const{
      return *(_allStraws_p.at(id.asUint16()));
    }

    std::array<Straw,StrawId::_nustraws> const& getStraws() const{
      return _allStraws;
    }

    bool strawExists( StrawId const id) const{
      // return _allStraws_p.at(id.asUint16()) != nullptr;
      return _strawExists2.at(id.asUint16());
    }

    // Return G4TUBS parameters for straws, includes
    // wire, gas and straw materials.
    TubsParams strawOuterTubsParams(StrawId const& id) const;
    TubsParams strawWallMother(StrawId const& id) const;
    TubsParams strawWallOuterMetal(StrawId const& id)  const;
    TubsParams strawWallInnerMetal1(StrawId const& id) const;
    TubsParams strawWallInnerMetal2(StrawId const& id) const;
    TubsParams strawWireMother(StrawId const& id) const;
    TubsParams strawWirePlate(StrawId const& id) const;

  protected:

    std::string _name;

    // Position of the center of the tracker, in the Mu2e coordinate system.
    double _z0;

    // Outer radius of a logical volume that will just contain the entire tracker.
    double _rOut;

    // All envelope volumes are made of this.
    std::string _envelopeMaterial;

    double _strawInnerRadius;
    double _strawOuterRadius;
    double _strawWallThickness;
    double _outerMetalThickness;
    double _innerMetal1Thickness;
    double _innerMetal2Thickness;
    double _wireRadius;
    double _wirePlateThickness;

    std::string _wallMaterialName;
    std::string _outerMetalMaterial;
    std::string _innerMetal1Material;
    std::string _innerMetal2Material;
    std::string _gasMaterialName;
    std::string _wireMaterialName;
    std::string _wirePlateMaterial;

    // Lengths of straws indexed by manifold, from innermost radius, outwards.
    std::array<double,StrawId::_nstraws> _strawHalfLengths;

    // Same for the active length of the straw.
    // This is only valid for SupportModel==simple
    std::array<double,StrawId::_nstraws> _strawActiveHalfLengths;

    // Deprecated: part of the ancient MECO Tracker design.
    // A few vestiges not yet removed.
    std::vector<Manifold> _allManifolds;

    // Outer envelope that holds the new style support structure.
    PlacedTubs _mother;

    // The envelope that holds all of the planes in the tracker,
    // including the plane supports.
    TubsParams _innerTrackerEnvelopeParams;

    // The envelope that holds all of the pieces in one plane, including supports.
    TubsParams _planeEnvelopeParams;

    // Ditto for Panel
    TubsParams _panelEnvelopeParams;

    // Which level of detail is present in the model of the support structure?
    SupportModel _supportModel;

    // All supports are the same shape; only relevant for _supportModel=="simple"
    Support _supportParams;
    double _panelZOffset; // introduced for version 5

    // A more detailed model of the supports; again each plane has identical supports.
    // only relevant for _supportModel == "detailedv0".
    SupportStructure _supportStructure;

    // All manifolds are the same shape.
    // Deprecated: these will go away soon.
    std::vector<double> _manifoldHalfLengths;

    // Inner radius of inside edge of innermost straw.
    double _envelopeInnerRadius;

    // presence info for each straw.
    //    std::vector<bool> _strawExists;

    // =============== NewTracker Private Objects Start ==============

    // Dense arrays
    std::array<Plane,StrawId::_nplanes> _planes;
    std::array<Panel,StrawId::_nupanels> _panels;
    std::array<Straw,StrawId::_nustraws> _allStraws;

    // Sparse array: designed for indexing by StrawId.
    // For all legal entries in StrawId, this points to a straw in _straws2;
    // All other entries are null.
    std::array<Straw const*,Tracker::_maxRedirect> _allStraws_p;

    // Another sparse array
    std::array<bool,Tracker::_maxRedirect> _strawExists2;

    // =============== NewTracker Private Objects End ==============


  private:
    // copying is complex
    // prevent these from being generated or used
    Tracker& operator=(const Tracker& other);
    Tracker(Tracker&& other) noexcept; // move constructor
    Tracker& operator=(Tracker&& other) noexcept; // move assignment

  };

} //namespace mu2e

#endif /* TrackerGeom_Tracker_hh */
