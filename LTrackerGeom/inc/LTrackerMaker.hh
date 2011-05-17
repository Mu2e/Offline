#ifndef LTrackerGeom_LTrackerMaker_hh
#define LTrackerGeom_LTrackerMaker_hh
//
// Construct and return an LTracker.
//
//
// $Id: LTrackerMaker.hh,v 1.6 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <memory>

#include "LTrackerGeom/inc/LayerInfo.hh"

#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/Sector.hh"
#include "TrackerGeom/inc/Layer.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LTracker;
  class SimpleConfig;

  class LTrackerMaker{

  public:
    LTrackerMaker( int nSides,
                   std::vector<LayerInfo> sideInfo,
                   std::vector<LayerInfo> vaneInfo,
                   double r0,
                   double halfLength,
                   double radius,
                   CLHEP::Hep3Vector center,
                   double phi0,
                   double tiltX,
                   double tiltY,
                   CLHEP::Hep3Vector vaneOffset
                   );

    LTrackerMaker( SimpleConfig const& config );  
    ~LTrackerMaker ();

    // This is depracted and will go away soon.  
    // Still needed for root graphics version.
    const LTracker& getLTracker() const { return *_ltt;}

    // This is the accessor that will remain.
    std::auto_ptr<LTracker> getLTrackerPtr() { return _ltt; }

  private:

    void BuildIt();
    void CheckFit();
    void CheckSideConsistency();
    void CheckVaneConsistency();
    void checkForOverlaps( bool printWarnings = true );
    int  totalStraws() const;

    void MakeSides();
    void MakeVanes();
    void MakeDetails();
    void FillNearestNeighbours();
    void FillPointersAndIndices();
    void FillPointersAndIndices2();

    //  SimpleConfig const& _config;

    // Number of sides in the polygon.
    int _nSides;

    // Number of straws per layer in sides.
    std::vector<LayerInfo> _sideInfo;

    // Number of straws per layer in vanes.
    std::vector<LayerInfo> _vaneInfo;

    // Nominal radius of the tracker; radius to mid-plane of a layer at z=0.
    double _r0;
  
    // Length of the tracker mother volume.
    double _halfLength;

    // Outerradius of a volume the tracker mother volume.
    double _rOut;

    // Passed through to LTracker as is.
    double _z0;
    std::string _fillMaterial;

    // Names of the materials in the straws.
    std::vector<std::string> _strawMaterialNames0, _strawMaterialNames1;

    // Straw half length, radius and thickness; thickness of carbon for conducting straws.
    double _strawHalfLength;
    double _strawRadius;
    double _strawThick;
    double _carbonThick;

    // Radius of the wire.
    double _rwire;

    // Center of the tracker.
    CLHEP::Hep3Vector _center;

    // Overall azimuthal rotation
    double _phi0;

    // See Mu2e-doc-561-v2 ( or higher version ).
    double _tiltX;
    double _tiltY;

    // 
    CLHEP::Hep3Vector _vaneOffset;

    // Half of the angle subtended by one side of the polygon.
    // Its tangent, cosine and sine.
    double _phiHalf;
    double _tphiHalf;
    double _cphiHalf;
    double _sphiHalf;

    std::auto_ptr<LTracker> _ltt;

  };

}  //namespace mu2e

#endif /* LTrackerGeom_LTrackerMaker_hh */
