#ifndef GeometryService_TTrackerMaker_hh
#define GeometryService_TTrackerMaker_hh
//
// Construct and return a TTracker.
//
// Original author Rob Kutschke
//

#include <memory>
#include <string>
#include <vector>

#include "TTrackerGeom/inc/TTracker.hh"

#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/Layer.hh"
#include "TrackerGeom/inc/Panel.hh"

#include "TTrackerGeom/inc/Station.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SimpleConfig;
  class TTracker;

  class TTrackerMaker{

  public:

    TTrackerMaker( SimpleConfig const& config );

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    std::unique_ptr<TTracker> getTTrackerPtr() { return std::move(_tt); }

  private:

    // Extract info from the config file.
    void parseConfig( const SimpleConfig& config );

    // Fill the details of the different straw types.
    void makeDetails();

    void makeMother();
    void makeDevice( DeviceId devId );
    void makePanel( const PanelId& secId, Device& dev );
    void makeLayer ( const LayerId& layId,  Panel& sec );
    void makeManifolds( const PanelId& secId);

    void computeStrawHalfLengths();
    void computePanelBoxParams(Panel& panel, Device& dev);
    void computeConstantPanelBoxParams();
    void computeLayerSpacingAndShift();
    void computeManifoldEdgeExcessSpace();
    void computeTrackerEnvelope();
    void computeDeviceEnvelope();

    void identifyNeighbourStraws();
    void identifyDirectionalNeighbourStraws();
    StrawIndex ttStrawIndex (LayerId const &layerId, int snum);

    void makeStation( StationId stationId );
    void makePlaneMF  ( const PlaneMFId &planeMFId, Station & station );
    void makeFace   ( const FaceId &faceId, PlaneMF &planeMF, const Device &device );
    void makePanelMF  ( const PanelMFId &panelMFId, Face &face, const Panel &panel );
    void makeZLayer ( const ZLayerId &zlayerId, PanelMF &panelMF, const Layer &layer );

    // Do the work of constructing it.
    void buildIt();

    double chooseDeviceSpacing( int idev ) const;
    double findFirstDevZ0() const;
    double panelRotation(int isec,int idev) const;

    // Some functions for the, fully detailed support structure.
    void makeSupportStructure();
    void makeStrawTubs();
    void makeThinSupportRings();

    // An ugly hack for the detailed support structure; must be called after all
    // straws have been created.
    void recomputeHalfLengths();

    // Some final self-consistency checks.
    void finalCheck();

    int    _verbosityLevel;

    // Properties of the mother volume.
    double _motherRIn;               // Inner radius of the mother volume
    double _motherROut;              // Outer radius of the mother volume
    double _motherHalfLength;        // Half length, in z, of the mother
    double _motherZ0;                // Mother volume is centered on the DS axis; the z of the center is NOT
                                     // necessarily the same as the center of the instrumented region of the tracker.

    // Basic parameters needed to describe the TTracker.
    int    _numDevices;                  // Number of devices.
    int    _panelsPerDevice;            // Number of panels in one device.
    int    _layersPerPanel;             // Number of layers in one panel.
    int    _manifoldsPerEnd;             // Number of manifolds along one end of the wires in a layer.
    int    _strawsPerManifold;           // Number of straws connected to each manifold.
    int    _rotationPattern;             // Pattern of rotations from device to device.
    int    _spacingPattern;              // Pattern of spacing from device to device.
    double _oddStationRotation;           // rotation of odd stations relative to even
    double _zCenter;                     // z position of the center of the tracker, in the Mu2e coord system.
    double _xCenter;                     // x position of the center of the tracker, in the Mu2e coord system.
    double _envelopeInnerRadius;         // Inner radius of inside of innermost straw.
    double _strawOuterRadius;            // Radius of each straw.
    double _strawWallThickness;          // Thickness of each straw.
    double _strawGap;                    // Gap between straws.
    double _deviceSpacing;               // Z-separation between adjacent stations.
    double _deviceHalfSeparation;        // Z-separation between adjacent devices.
    double _innerSupportRadius;          // Inner radius of support frame.
    double _outerSupportRadius;          // Outer radius of support frame.
    double _supportHalfThickness;        // Thickness of support frame.
    double _wireRadius;                  // Wire radius
    double _manifoldYOffset;             // Offset of manifold from inner edge of support.
    double _passivationMargin;           // Deadened region near the end of the straw.
    double _wallOuterMetalThickness;     // Metal on the outer surface of the straw
    double _wallInnerMetal1Thickness;    // Metal on the inner surface of the straw
    double _wallInnerMetal2Thickness;    // Thinner metal layer between Metal1 and the gas
    double _wirePlateThickness;          // Plating on the top of the wire.
    double _virtualDetectorHalfLength;   // Half length (or thickness) of virtual detectors - needed to leave space for them!

    std::string _envelopeMaterial;        // Material for the envelope volume.
    std::string _supportMaterial;         // Material for the support volume.
    std::string _wallOuterMetalMaterial;  // Material on the outer surface of the straw
    std::string _wallInnerMetal1Material; // Material on the inner surface of the straw
    std::string _wallInnerMetal2Material; // Material between Meta1 and the gas.
    std::string _wirePlateMaterial;       // Material for plating on the wire.

    std::vector<double> _manifoldHalfLengths; // Dimensions of each manifold.
    std::vector<std::string> _strawMaterials; // Names of the materials.

    // Base rotations of a panel; does not include device rotation.
    std::vector<double> _panelBaseRotations;
    std::vector<double> _panelZSide;
    double _devrot; // hack to make redundant information self-consistent

    std::unique_ptr<TTracker> _tt;

    // Derived parameters.

    // Compute at start and compare against final result.
    int _nStrawsToReserve;

    // Lengths of straws indexed by manifold, from innermost radius, outwards.
    std::vector<double> _strawHalfLengths;

    // Same for the active length of the straw.
    // This is only valid for SupportModel==simple
    std::vector<double> _strawActiveHalfLengths;

    // panel box half lengths
    std::vector<double> _panelBoxHalfLengths;

    // panel box offset magnitudes
    double _panelBoxXOffsetMag;
    double _panelBoxZOffsetMag;

    // distance between layers (in Z)
    double _layerHalfSpacing;
    // relative layer shift (in X)
    double _layerHalfShift;

    // Space between first/last straw and edge of manifold
    double _manifoldXEdgeExcessSpace;
    double _manifoldZEdgeExcessSpace;

    // Z Location of the first device.
    double _z0;

    // Parameters used in the Station/PlaneMF/Face/PanelMF/ZLayer view
    int    _numStations;                  // Number of Stations.
    int    _planeMFsPerStation;
    int    _facesPerStation;
    int    _facesPerPlaneMF;
    int    _panelMFsPerFace;
    int    _zLayersPerPanelMF;
    int    _strawsPerZLayer;

    // The detailed description of the complete support structure
    SupportModel _supportModel;
    double       _endRingOuterRadius;
    double       _endRingInnerRadius;
    double       _endRingHalfLength;
    double       _endRingZOffset;
    std::string  _endRingMaterial;
    std::vector<int> _midRingSlot;
    double           _midRingHalfLength;

    double      _innerRingInnerRadius;
    double      _innerRingOuterRadius;
    double      _innerRingHalfLength;
    std::string _innerRingMaterial;

    double      _centerPlateHalfLength;
    std::string _centerPlateMaterial;

    double      _outerRingInnerRadius;
    double      _outerRingOuterRadius;
    std::string _outerRingMaterial;

    double      _coverHalfLength;
    std::string _coverMaterial;

    double      _electronicsG10HalfLength;
    std::string _electronicsG10Material;

    double      _electronicsCuHhalfLength;
    std::string _electronicsCuMaterial;

    double      _channelZOffset;
    double      _channelDepth;
    std::string _channelMaterial;

    std::string _electronicsSpaceMaterial;

    int         _nStaves;
    double      _stavePhi0;
    double      _stavePhiWidth;
    std::string _staveMaterial;

    std::vector<double> _beam0_phiRange;
    double _beam0_innerRadius;
    double _beam0_outerRadius;
    std::string _beam0_material;

    std::vector<double> _beam1_phiRange;
    std::vector<double> _beam1_phiSpans;
    std::vector<double> _beam1_servicePhi0s;
    std::vector<double> _beam1_servicePhiEnds;
    double _beam1_innerRadius;
    double _beam1_midRadius1;
    std::vector<double> _beam1_serviceOuterRadii;
    double _beam1_midRadius2;
    double _beam1_outerRadius;
    std::string _beam1_material;
    std::vector<std::string> _beam1_serviceMaterials;
    std::vector<double> _beam1_serviceCovRelThickness;
    std::vector<std::string> _beam1_serviceMaterialsCov;

    std::vector<int>  _nonExistingDevices;

  };

}  //namespace mu2e

#endif /* GeometryService_TTrackerMaker_hh */
