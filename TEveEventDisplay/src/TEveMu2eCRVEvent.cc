#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
using namespace mu2e;
namespace mu2e{

	TEveMu2eCRVEvent::TEveMu2eCRVEvent(){}
  
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
	{
    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(mColor);
    this->SetMarkerSize(mSize);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }

  void TEveMu2eCRVEvent::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(mColor);
    this->SetMarkerSize(mSize);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }
}

