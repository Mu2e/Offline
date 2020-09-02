#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCTraj.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
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

 void TEveMu2eMCTraj::DrawLine3D(const std::string &pstr,  CLHEP::Hep3Vector Start, CLHEP::Hep3Vector End, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, -1)).c_str());
    //hep3vectorTocm(Start);
    std::cout<<" Start "<<Start.x()<<" "<<Start.y()<<" "<<Start.z()<<std::endl;
    //hep3vectorTocm(End);
    std::cout<<" End "<<End.x()<<" "<<End.y()<<" "<<End.z()<<std::endl;
    TEveLine *line = new TEveLine();
    line->SetPoint(0, Start.x()/10-390.4, Start.y()/10, Start.z()/10+1017.5-TrackerLength()/2); //y+231.2
    line->SetNextPoint(End.x()/10-390.4, End.y()/10, End.z()/10+1017.5-TrackerLength()/2); //1017.5 - halflength
    line->SetLineColor(kYellow);
    this->SetMarkerColor(kBlue);
    this->SetMarkerSize(2);
    this->SetPickable(kTRUE);
    HitList->AddElement(line);
    HitList->AddElement(this);
  }
}
