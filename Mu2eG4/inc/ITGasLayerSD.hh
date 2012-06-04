#ifndef Mu2eG4_ITGasLayerSD_hh
#define Mu2eG4_ITGasLayerSD_hh
//
// Defines a generic ITracker sensitive detector
//
// $Id: ITGasLayerSD.hh,v 1.10 2012/06/04 23:46:23 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/06/04 23:46:23 $
//
// Original author G. Tassielli
//

//// Mu2e includes
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  class ITGasLayerSD : public Mu2eSensitiveDetector {

  public:
    ITGasLayerSD(G4String, const SimpleConfig& config);
    ~ITGasLayerSD();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) { return false; }

    static void setMu2eDetCenterInWorld(const G4ThreeVector &origin) {
            _mu2eDetCenter = origin;
    }

  protected:

    int _superlayer;
    int _ring;
    int _nwires;
    double _Dphi;

    GeomHandle<ITracker> itracker;
    //ITracker::GeomType _ittype;

    // Mu2e point of origin
    static G4ThreeVector _mu2eDetCenter;

  };

} // namespace mu2e

#endif /* Mu2eG4_ITGasLayerSD_hh */
