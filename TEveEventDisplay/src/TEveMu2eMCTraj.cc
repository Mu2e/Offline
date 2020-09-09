#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCTraj.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
using namespace mu2e;
namespace mu2e{

  void TEveMu2eMCTraj::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(kBlue);
    this->SetMarkerSize(2);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }

 void TEveMu2eMCTraj::DrawLine(const std::string &pstr,  CLHEP::Hep3Vector Start, CLHEP::Hep3Vector End, TEveElementList *HitList)
  {
    
    this->SetTitle((DataTitle(pstr, -1)).c_str());
   
    hep3vectorTocm(Start);
    hep3vectorTocm(End);
    TEveLine *line = new TEveLine();
 
    line->SetPoint(0, Start.x(), Start.y(), Start.z()); 
    line->SetNextPoint(End.x(), End.y(), End.z()); 
    
    line->SetLineColor(kYellow);
    this->SetMarkerColor(kYellow);
    this->SetMarkerSize(5);
    this->SetPickable(kTRUE);
    HitList->AddElement(line);
    HitList->AddElement(this);
  }
}
