#ifndef TrackerGeom_TrackerG4Info_hh
#define TrackerGeom_TrackerG4Info_hh
///
// Tracker geometry content specific to building the G4 model
// Extracted from the original Tracker
//
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/GeomPrimitives/inc/PlacedTubs.hh"
#include "Offline/TrackerGeom/inc/PanelEB.hh"
#include "Offline/TrackerGeom/inc/Manifold.hh"
#include "Offline/TrackerGeom/inc/Support.hh"
#include "Offline/TrackerGeom/inc/SupportModel.hh"
#include "Offline/TrackerGeom/inc/SupportStructure.hh"

namespace mu2e {
  class TrackerG4Info {
    friend class TrackerMaker;
    public:
    // electronics board
    const PanelEB& panelElectronicsBoard() const { return _panelEB;}
    const SupportModel& getSupportModel() const{ return _supportModel; }
    const Support& getSupportParams () const{ return _supportParams; }
    const SupportStructure& getSupportStructure() const{ return _supportStructure; }
    const TubsParams& getPlaneEnvelopeParams() const{ return _planeEnvelopeParams; }
    const TubsParams& getPanelEnvelopeParams() const{ return _panelEnvelopeParams; }
    const TubsParams& getInnerTrackerEnvelopeParams() const{ return _innerTrackerEnvelopeParams; }
    const PlacedTubs& mother() const{ return _mother; }
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

    double panelOffset() const { return _panelZOffset; }
    double z0()   const { return _z0;} // in Mu2e coordinates

    private:
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
    // Electronics board
    PanelEB _panelEB;

    double _panelZOffset; // introduced for version 5
    // Position of the center of the tracker, in the Mu2e coordinate system.
    double _z0;

  };
}
#endif
