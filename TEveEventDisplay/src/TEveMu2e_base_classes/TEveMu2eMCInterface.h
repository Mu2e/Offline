#ifndef TEveMu2eMCInterface_h
#define TEveMu2eMCInterface_h

//libGeom
#include <TGeoManager.h>
//TEve
#include <TEveManager.h>
#include <TEveStraightLineSet.h>
//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "MCDataProducts/inc/MCTrajectoryPoint.hh"
//TEveMu2e
#include "TEveEventDisplay/src/dict_classes/Collection_Filler.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCTraj.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
namespace mu2e{
    class TEveMu2eMCInterface {
    public:
      #ifndef __CINT__
      TEveMu2eMCInterface() : fTrackList2D(0),fTrackList3D(0){};
      TEveMu2eMCInterface(const TEveMu2eMCInterface &);
      TEveMu2eMCInterface& operator=(const TEveMu2eMCInterface &);
      virtual ~TEveMu2eMCInterface(){};
      void AddSimpleMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);
      void AddFullMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);
      #endif
      TEveElementList *fTrackList2D;
      TEveElementList *fTrackList3D;

      ClassDef(TEveMu2eMCInterface,0);

  }; //end class def

}//end namespace mu2e

#endif /*TEveMu2eMCInterface.h*/
