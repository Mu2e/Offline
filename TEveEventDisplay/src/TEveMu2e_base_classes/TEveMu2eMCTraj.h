#ifndef TEveMu2eMCTraj_h
#define TEveMu2eMCTraj_h

#include <TObject.h>
#include <string.h>
#include <TEvePointSet.h>
#include <TEveLine.h>
#include <TPolyLine3D.h>
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
namespace mu2e {
  class TEveMu2eMCTraj : public TEvePointSet{
    public:
      #ifndef __CINT__
      explicit TEveMu2eMCTraj();
      virtual ~TEveMu2eMCTraj(){};
      void DrawHit3D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list);
      void DrawSimpleLine(const std::string &pstr, CLHEP::Hep3Vector Start, CLHEP::Hep3Vector End, TEveElementList *HitList);
      void DrawFullLine(const std::string &pstr, CLHEP::Hep3Vector Start, CLHEP::Hep3Vector End, TEveElementList *HitList);
      std::string DataTitle(const std::string &pstr, Int_t n);
      #endif
      ClassDef(TEveMu2eMCTraj, 0);
    };
}
#endif
