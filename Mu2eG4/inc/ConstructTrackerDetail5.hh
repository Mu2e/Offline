#ifndef Mu2eG4_ConstructTrackerDetail5_hh
#define Mu2eG4_ConstructTrackerDetail5_hh
//
// Class to construct the detailed version of the Tracker as of Nov 2017
//
//
// Based on Original code by Rob Kutschke and Krzysztof Genser
// This implementation by David N. Brown (Louisville), November 2017
//

#include "G4Helper/inc/VolumeInfo.hh"

#include "G4RotationMatrix.hh"

#include <vector>

namespace mu2e {

  class SimpleConfig;
  class G4Helper;
  class AntiLeakRegistry;
  class Tracker;

  class ConstructTrackerDetail5 {

  public:
    ConstructTrackerDetail5( VolumeInfo   const& ds3Vac,
                              SimpleConfig const& config );

    VolumeInfo motherInfo() { return _motherInfo; }


  private:

    CLHEP::Hep3Vector computeOffset();

    // Mother volume to hold all other G4 geom objects created by this class.
    void constructMother();

    // Construct the main support structure: initially this is the
    // stiffening rings and staves.
    void constructMainSupports();

    // Construct all Planes.
    void constructPlanes();

    // Construct all panels(panels) within one Plane
    void addPanelsAndEBKeys( VolumeInfo&     baseStrawPanel,
			     int              idev,
                             VolumeInfo&      plane,
                             VolumeInfo&      baseEBKey,
                             VolumeInfo&      baseEBKeyShield,
                             VolumeInfo&      trackerMother,
                             double           panelCenterPhi );


    // Build logical volume heirarchy for one panel: straws placed inside 
    // a panel mother volume with supports, then placed in a plane
    VolumeInfo preparePanel(const int& ipln, const int& ipnl,
			    VolumeInfo& thePlane, VolumeInfo& strawPanel,
			    CLHEP::Hep3Vector& pnlPos,
			    G4RotationMatrix* rot);

    // This is like the former preparePanel function - it just makes a set of
    // straws for repeated use in panels.
    VolumeInfo prepareStrawPanel();

    VolumeInfo prepareEBKey(bool keyItself);

    inline double  getWithinZeroTwoPi (double phi0);
    inline double diffWithinZeroTwoPi (double phi1, double phi2);

    // Compute the half width, in phi, of a panel envelope.
    double panelHalfAzimuth();

    // For debubbing, construct boxes that represent coordinate axes for debugging.
    void constructAxes();

    // References to c'tor arguments.
    VolumeInfo   const& _ds3Vac;
    SimpleConfig const& _config;

    G4Helper&          _helper;
    AntiLeakRegistry & _reg;

    Tracker const& _tracker;

    int  _verbosityLevel;
    bool _doSurfaceCheck;
    bool _forceAuxEdgeVisible;

    // Position of tracker center inside of its mother volume.
    CLHEP::Hep3Vector _offset;

    // Information about the mother volume.
    VolumeInfo _motherInfo;

  }; // end class ConstructTrackerDetail5

} // end namespace mu2e

#endif
