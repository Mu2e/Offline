#ifndef TEveMu2eTracker_h
#define TEveMu2eTracker_h

#include <TApplication.h>
#include<TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <TEveGeoShape.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
// ... libRIO
#include <TFile.h>
//Tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
//TEveMu2e:
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"

namespace mu2e {
	
class TEveMu2eTracker{
	public:
    #ifndef __CINT__
    explicit TEveMu2eTracker();
    TEveMu2eTracker(const TEveMu2eTracker &){};
    TEveMu2eTracker& operator=(const TEveMu2eTracker &);
    virtual ~TEveMu2eTracker(){};
    TEveMu2e2DProjection *tracker2Dproj = new TEveMu2e2DProjection();
    void DrawTrackerDetector(art::Run const& run, TGeoVolume* topvol , TEveElementList *orthodet);
    #endif
    ClassDef(TEveMu2eTracker, 0);
  };
}
#endif
