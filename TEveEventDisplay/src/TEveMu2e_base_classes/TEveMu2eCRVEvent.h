#ifndef TEveMu2eCRVEvent_h
#define TEveMu2eCRVEvent_h

#include <TObject.h>

#include <TEvePointSet.h>
#include <TEveLine.h>
#include "RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
 
namespace mu2e {
  class TEveMu2eCRVEvent : public TEvePointSet {
    public:
      #ifndef __CINT__
      explicit TEveMu2eCRVEvent();
      TEveMu2eCRVEvent(CrvRecoPulse chit) : fCrvRecoPulse(chit){};
      virtual ~TEveMu2eCRVEvent(){};

      CrvCoincidenceCluster fCrvCoincidenceCluster; 
      CrvRecoPulse fCrvRecoPulse;

      Int_t mColor = kBlue;
      Int_t mSize= 1; 

      bool AddErrorBar = true;
      void DrawHit2D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list); 
      void DrawHit3D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list); 
      std::string DataTitle(const std::string &pstr, int n);
      #endif
      ClassDef(TEveMu2eCRVEvent, 0);
  };
  typedef std::vector<mu2e::TEveMu2eCRVEvent> TEveMu2eCRVEventCollection;
}
#endif
