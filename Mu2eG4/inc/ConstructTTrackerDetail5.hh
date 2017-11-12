#ifndef Mu2eG4_ConstructTTrackerDetail5_hh
#define Mu2eG4_ConstructTTrackerDetail5_hh
//
// Class to construct the detailed version of the TTracker as of Nov 2017
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
  class TTracker;

  class ConstructTTrackerDetail5 {

  public:
    ConstructTTrackerDetail5( VolumeInfo   const& ds3Vac,
                          SimpleConfig const& config  );

    VolumeInfo motherInfo() { return _motherInfo; }


  private:

    CLHEP::Hep3Vector computeOffset();

    // Mother volume to hold all other G4 geom objects created by this class.
    void constructMother();

    // Construct the main support structure: initially this is the
    // stiffening rings and staves.
    void constructMainSupports();

    // Construct all stations.
    void constructStations();

    // Construct all panels(panels) within one station
    void addPanelsAndEBKeys( int              idev,
                             VolumeInfo&      plane,
                             VolumeInfo&      baseEBKey,
                             VolumeInfo&      baseEBKeyShield,
                             VolumeInfo&      trackerMother,
                             double           panelCenterPhi );

    // Construct the support infrastructure for each station.
    void addPlaneSupports( std::vector<VolumeInfo>& supportsInfo, int idev, VolumeInfo const& devInfo );

    // Build logical volume heirarchy for the elements of the support structure that are inside
    // each plane(plane) envelope.  Do not place this volume heirarchy.
    //    void preparePlaneSupports( std::vector<VolumeInfo>& supportsInfo );

    // Build logical volume heirarchy for one panel: straws placed inside a panel mother volume.
    // Do not place this volume heirarchy.
    VolumeInfo preparePanel(const int& ipln, const int& ipnl,
			    VolumeInfo& thePlane, CLHEP::Hep3Vector& pnlPos,
			    G4RotationMatrix* rot);

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

    TTracker const& _ttracker;

    int  _verbosityLevel;
    bool _doSurfaceCheck;
    bool _forceAuxEdgeVisible;

    // Position of tracker center inside of its mother volume.
    CLHEP::Hep3Vector _offset;

    // Information about the mother volume.
    VolumeInfo _motherInfo;

  }; // end class ConstructTTrackerDetail5

} // end namespace mu2e

#endif
