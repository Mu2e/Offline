#ifndef TEveMu2eCRVEvent_h
#define TEveMu2eCRVEvent_h

#include <TObject.h>

#include <TEvePointSet.h>
#include <TEveLine.h>
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"

namespace mu2e {
  class TEveMu2eCRVEvent : public TEvePointSet {
    public:
      #ifndef __CINT__
      explicit TEveMu2eCRVEvent();
      TEveMu2eCRVEvent(CrvRecoPulse chit) : fCrvRecoPulse_(chit){};
      virtual ~TEveMu2eCRVEvent(){};

      CrvCoincidenceCluster fCrvCoincidenceCluster_;
      CrvRecoPulse fCrvRecoPulse_;

      Int_t mColor_ = kBlue;
      Int_t mSize_ = 1;
      bool AddErrorBar_ = true;
      std::tuple<CLHEP::Hep3Vector, CLHEP::Hep3Vector, std::string, int> DrawSciBar();
      void DrawHit2D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list2DXY, TEveElementList *list2DXZ);
      void DrawHit3D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list3D);
      std::string DataTitle(const std::string &pstr, int n);
      #endif
      ClassDef(TEveMu2eCRVEvent, 0);
  };
  typedef std::vector<mu2e::TEveMu2eCRVEvent> TEveMu2eCRVEventCollection;
}
#endif


