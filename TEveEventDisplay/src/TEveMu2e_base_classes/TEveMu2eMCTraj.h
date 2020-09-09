#ifndef TEveMu2eMCTraj_h
#define TEveMu2eMCTraj_h

#include <TObject.h>
#include <string.h>
#include <TEvePointSet.h>
#include <TEveLine.h>
#include <TPolyLine3D.h>
#include "RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {
  class TEveMu2eMCTraj : public TEvePointSet{//, public TPolyLine3D {
    public:
      #ifndef __CINT__
      explicit TEveMu2eMCTraj(){};
      virtual ~TEveMu2eMCTraj(){};
      void DrawHit3D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector HitPos, TEveElementList *list); 
      void DrawLine(const std::string &pstr, CLHEP::Hep3Vector Start, CLHEP::Hep3Vector End, TEveElementList *HitList);
      inline std::string DataTitle(const std::string &pstr, Int_t n){
        std::string dstr = "";
        if (n != -1){dstr=" hit#" + std::to_string(n) + "\nLayer: ";}
        std::string strlab=pstr+dstr;
        return (strlab);
      }
      #endif
      ClassDef(TEveMu2eMCTraj, 0);
    };
}
#endif
