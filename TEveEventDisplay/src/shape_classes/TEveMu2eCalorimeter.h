#ifndef TEveMu2eCalorimeter_h
#define TEveMu2eCalorimeter_h

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
//Calorimeter:
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "CalorimeterGeom/inc/CaloInfo.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskGeomInfo.hh"
#include "CalorimeterGeom/inc/Disk.hh"
//TEveMu2e:
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"

namespace mu2e{
	class TEveMu2eCalorimeter
	{
    public:
      #ifndef __CINT__
      explicit TEveMu2eCalorimeter();
      TEveMu2eCalorimeter(const TEveMu2eCalorimeter &){};
      TEveMu2eCalorimeter& operator=(const TEveMu2eCalorimeter &);
      virtual ~TEveMu2eCalorimeter(){};
      void DrawCaloDetector(art::Run const& run, TGeoVolume* topvol , TEveElementList *orthodet0,TEveElementList *orthodet1);
      TEveMu2e2DProjection *calo2Dproj = new TEveMu2e2DProjection();
      #endif
      ClassDef(TEveMu2eCalorimeter, 0);
  };
}
#endif
