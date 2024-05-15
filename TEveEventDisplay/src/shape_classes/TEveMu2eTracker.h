#ifndef TEveMu2eTracker_h
#define TEveMu2eTracker_h

#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <TEveGeoShape.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
#include <TBox.h>
#include <TGeoBBox.h>
// ... libRIO
#include <TFile.h>
//Tracker
#include "Offline/TrackerGeom/inc/Tracker.hh"
//TEveMu2e:
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "Offline/TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/TargetFoil.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include <string>

namespace mu2e {

class TEveMu2eTracker{
  public:
    #ifndef __CINT__
    explicit TEveMu2eTracker();
    TEveMu2eTracker(const TEveMu2eTracker &){};
    TEveMu2eTracker& operator=(const TEveMu2eTracker &);
    virtual ~TEveMu2eTracker(){};
    TEveMu2e2DProjection *tracker2Dproj = new TEveMu2e2DProjection();
    void DrawTrackerDetector(TGeoVolume* topvol , TEveElementList *orthodetXZ, TEveElementList *orthodetXY);
    #endif
  ClassDef(TEveMu2eTracker, 0);
  };
}
#endif
