#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
#include "Offline/TEveEventDisplay/src/dict_classes/GeomUtils.h"
using namespace mu2e;
namespace mu2e{

	TEveMu2eCRVEvent::TEveMu2eCRVEvent(){}
	
	/*------------Function to build title:-------------*/
  std::string TEveMu2eCRVEvent::DataTitle(const std::string &pstr, int n){
    std::string dstr=" hit#" + std::to_string(n) + "\nLayer: ";
    std::string strlab=pstr+dstr;
    return (strlab);
  }

  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }
  
  /*------------Function to 2D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }
}

