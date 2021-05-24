#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"

using namespace mu2e;
namespace mu2e{

  TEveMu2eHit::TEveMu2eHit(){}
  
  /*------------Function to build title:-------------*/
  std::string TEveMu2eHit::DataTitle(const std::string &pstr, int n){
        std::string dstr=" hit#" + std::to_string(n) + "\nLayer: ";
        std::string strlab=pstr+dstr;
        return (strlab);
  }

  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eHit::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    if(AddErrorBar_){ 
      TEveLine *error = new TEveLine();
      auto const& p = fComboHit_.pos();
      auto const& w = fComboHit_.wdir();
      auto const& s = fComboHit_.wireRes();
      double x1 = (p.x()+s*w.x());
      double x2 = (p.x()-s*w.x());
      double z1 = (p.z()+s*w.z());
      double z2 = (p.z()-s*w.z());
      double y1 = (p.y()+s*w.y());
      double y2 = (p.y()-s*w.y());
      std::string errorbar = "ErrorBar Length: %d, %d, %d"; 
      error->SetTitle(Form(errorbar.c_str(), (x1 - x2), (y1 - y2), (z1 - z2)));
      GeomHandle<DetectorSystem> det;
      Hep3Vector vec1(x1, y1, z1);
      Hep3Vector vec2(x2, y2, z2);
      Hep3Vector inDet1 = det->toMu2e(vec1);
      Hep3Vector inDet2 = det->toMu2e(vec2);
      error->SetPoint(0, pointmmTocm(inDet1.x()),pointmmTocm(inDet1.y()),pointmmTocm(inDet1.z()));
      error->SetNextPoint(pointmmTocm(inDet2.x()), pointmmTocm(inDet2.y()),pointmmTocm(inDet2.z()));
      error->SetLineColor(kSpring);
      error->SetPickable(kTRUE);
      HitList->AddElement(error);
    }
    HitList->AddElement(this);
  }

  /*------------Function to 2D draw hits:-------------*/
  void TEveMu2eHit::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList)
  {
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);

    if(AddErrorBar_){ 
      TEveLine *error = new TEveLine();
      auto const& p = fComboHit_.pos();
      auto const& w = fComboHit_.wdir();
      auto const& s = fComboHit_.wireRes();
      double x1 = (p.x()+s*w.x());
      double x2 = (p.x()-s*w.x());
      double z1 = (p.z()+s*w.z());
      double z2 = (p.z()-s*w.z());
      double y1 = (p.y()+s*w.y());
      double y2 = (p.y()-s*w.y());
      
      std::string errorbar = "ErrorBar Length: %d, %d, %d"; 
      error->SetTitle(Form(errorbar.c_str(), (x1 - x2), (y1 - y2), (z1 - z2)));
      error->SetPoint(0, pointmmTocm(x1),pointmmTocm(y1),pointmmTocm(z1));
      error->SetNextPoint(pointmmTocm(x2),pointmmTocm(y2),pointmmTocm(z2));
      error->SetLineColor(kRed);
      error->SetPickable(kTRUE);
      HitList->AddElement(error);
    }
    HitList->AddElement(this);
    }
  }
