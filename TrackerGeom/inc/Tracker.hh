#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Hold all geometry and identifier information about
// a Tracker.  This is intended as a "data only" class.
//  Note that the only 'aligned' information is implicit in the individual straws
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
#include <memory>

#include "cetlib_except/exception.h"

#include "Mu2eInterfaces/inc/Detector.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "TrackerGeom/inc/Manifold.hh"
#include "TrackerGeom/inc/Support.hh"
#include "TrackerGeom/inc/SupportModel.hh"
#include "TrackerGeom/inc/SupportStructure.hh"
#include "GeneralUtilities/inc/HepTransform.hh"

#include "TrackerGeom/inc/Plane.hh"
#include "TrackerGeom/inc/Panel.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PlacedTubs.hh"
#include "TrackerGeom/inc/PanelEB.hh"
#include "TrackerGeom/inc/StrawProperties.hh"

namespace mu2e {
  class Tracker : public Detector, public ProditionsEntity {
    friend class TrackerMaker; // remove after factorizing out the G4 stuff FIXME!
    using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

    public:
    typedef std::shared_ptr<Tracker> ptr_t;
    typedef std::shared_ptr<const Tracker> cptr_t;
    using PlaneCollection = std::array<Plane,StrawId::_nplanes>;
    using PanelCollection = std::array<Panel,StrawId::_nupanels>;
    using StrawCollection = std::array<Straw,StrawId::_nustraws>;
    using StrawIndexMap = std::array<uint16_t,StrawId::_maxval>; 
    // default constructor results in non-functional object, but is required by proditions service
    Tracker() : _name("AlignedTracker") {}
    // copy constructor
    Tracker(const Tracker& other);
    // construct from a set of straws and their global properties.  It would be better to take ownership of
    // existing straws via std::move, but the need for the copy constructor above precludes that
    Tracker(StrawCollection const& straws, StrawProperties const& sprops, const char* name);
    // accessors

    std::string const& name() const { return _name; }

    // geometry parameters used only by G4.  These should be factorized out FIXME!
    double z0()   const { return _z0;} // in Mu2e coordinates (everything else is in Tracker coordiates).
    double zHalfLength() const;

    // Tracker Element Accessors
    size_t nPlanes() const { return _planes.size(); }
    size_t nPanels() const { return _panels.size(); }
    size_t nStraws() const { return _straws.size(); }

    // origin: defines nominal tracker coordinate system
    const xyzVec& origin() const { return _origin; }

    // constituent access
    const Plane& plane( const StrawId& id ) const{
      return _planes.at(id.getPlane());
    }

    PlaneCollection const& planes() const {
      return _planes;
    }
    PanelCollection const& panels() const{
      return _panels;
    }
    StrawCollection const& straws() const{
      return _straws;
    }

    const Panel& panel( const StrawId& id ) const{
      return _panels.at(id.uniquePanel());
    }

    uint16_t strawIndex(const StrawId& id) const {
      return _strawindex.at(id.asUint16());
    }

    const Straw& straw( const StrawId& id) const{
      return _straws[strawIndex(id)];
    }

// deprecated interface 
    const Plane& getPlane( const StrawId& id ) const{
      return _planes.at(id.getPlane());
    }
    
    const Plane& getPlane( uint16_t n ) const{
      return _planes.at(n);
    }

    PlaneCollection const& getPlanes() const {
      return _planes;
    }

    const Panel& getPanel( const StrawId& id ) const{
      return _panels.at(id.uniquePanel());
    }

    const Straw& getStraw( const StrawId& id) const{
      return _straws[strawIndex(id)];
    }

    Straw& getStraw( const StrawId& id) {
      return _straws[strawIndex(id)];
    }
    
    StrawCollection const& getStraws() const{
      return _straws;
    }


    // common attributes that don't depend on the specific element
    const StrawProperties& strawProperties() const { return _strawprops; }
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
    bool planeExists(StrawId const& id) const { return _planeExists[id.plane()]; }

    // G4 stuff: this should be factorized out TODO!
    std::string const& wallMaterialName()            const{ return _wallMaterialName; }
    std::string const& wallCoreMaterialName()        const{ return  wallMaterialName();  }
    std::string const& wallOuterMetalMaterialName()  const{ return _outerMetalMaterial;  }
    std::string const& wallInnerMetal1MaterialName() const{ return _innerMetal1Material; }
    std::string const& wallInnerMetal2MaterialName() const{ return _innerMetal2Material; }
    std::string const& gasMaterialName()             const{ return _gasMaterialName; }
    std::string const& wireMaterialName()            const{ return _wireMaterialName; }
    std::string const& wireCoreMaterialName()        const{ return  wireMaterialName();  }
    std::string const& wirePlateMaterialName()       const{ return _wirePlateMaterial;   }
    std::string const& envelopeMaterial()	     const { return _envelopeMaterial; }

//Mu2eG4 specific interface: these should be factored out TODO

    // electronics board
    const PanelEB& panelElectronicsBoard() const { return _panelEB;}
//    double rOut() const { return _rOut;}
    double panelOffset() const { return _panelZOffset; }
// why does this return by value??? FIXME!
    SupportModel getSupportModel() const{
      return _supportModel;
    }

    const Support& getSupportParams () const{
      return _supportParams;
    }

    const SupportStructure& getSupportStructure() const{
      return _supportStructure;
    }
// why do these return by value????? FIXME
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
  // G4-specific variables, set in TrackerMaker.  These should all go to G4Tracker FIXME!
    // Position of the center of the tracker, in the Mu2e coordinate system.
    double _z0;
    // Outer radius of a logical volume that will just contain the entire tracker.
    double _rOut; // is this ever used?  I think not

    // these should be in a struct FIXME!
    std::string _wallMaterialName;
    std::string _outerMetalMaterial;
    std::string _innerMetal1Material;
    std::string _innerMetal2Material;
    std::string _gasMaterialName;
    std::string _wireMaterialName;
    std::string _wirePlateMaterial;
    std::string _envelopeMaterial;

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

    // A more detailed model of the supports; again each plane has identical supports.
    // only relevant for _supportModel == "detailedv0".
    SupportStructure _supportStructure;

    // Deprecated: these will go away soon.
    std::vector<double> _manifoldHalfLengths;
    double _panelZOffset; // introduced for version 5

    // Inner radius of inside edge of innermost straw.
    double _envelopeInnerRadius;

   // Electronics board
    PanelEB _panelEB;

    // non-G4 content; this is the core of the class.  These are kept when the G4 content is factorized away to a separate class FIXME!
  private:
    std::string _name; // needed for proditions
    xyzVec _origin;
    HepTransform _tracker_to_tnom; // transform to NOMINAL tracker coordinate system.  This is a null transform for the nominal tracker
    // global straw properties
    StrawProperties _strawprops;
    // Dense arrays
    PlaneCollection _planes;
    PanelCollection _panels;
    // indirection from StrawId, for efficient lookup
    StrawIndexMap _strawindex;
    // plane existence: use cases of this should switch to using TrackerStatus and this should be removed FIXME!!
    std::array<bool,StrawId::_nplanes> _planeExists;
    // fundamental content is in the following
    StrawCollection _straws;
  };

} //namespace mu2e

#endif /* TrackerGeom_Tracker_hh */
