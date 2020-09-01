#ifndef TEveMu2eCRV_h
#define TEveMu2eCRV_h

#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <TEveGeoShape.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
// ... libRIO
#include <TFile.h>
//CRV/CRS:
#include "GeometryService/inc/GeomHandle.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
//TEveMu2e:
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"

namespace mu2e{
	class TEveMu2eCRV 
	{
    public:
      #ifndef __CINT__
      explicit TEveMu2eCRV();
      TEveMu2eCRV(const TEveMu2eCRV &){};
      TEveMu2eCRV& operator=(const TEveMu2eCRV &);
      virtual ~TEveMu2eCRV(){};
      void DrawCRVDetector(art::Run const& run, TGeoVolume* topvol , TEveElementList *orthodetlist[]);
      TEveMu2e2DProjection *CRV2Dproj = new TEveMu2e2DProjection();
      #endif
      ClassDef(TEveMu2eCRV, 0);
  };
}
#endif
