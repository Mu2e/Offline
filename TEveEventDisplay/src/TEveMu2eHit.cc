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
    //hep3vectorTocm(pointInMu2e);
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
  void TEveMu2eHit::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList2DXY, TEveElementList *HitList2DXZ)
  { 
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    const auto& allStraws = tracker->getStraws();
    for (size_t i = 0; i<tracker->nStraws(); i++)
      {
      const mu2e::Straw& s = allStraws[i];
      const CLHEP::Hep3Vector& p = s.getMidPoint();
      const CLHEP::Hep3Vector& d = s.getDirection();
      double x = p.x();
      double y = p.y();
      double z = p.z();
      double theta = d.theta();
      double phi = d.phi();
      double l = s.halfLength();
      double st=sin(theta);
      double ct=cos(theta);
      double sp=sin(phi+TMath::Pi()/2.0);   
      double cp=cos(phi+TMath::Pi()/2.0);
      double x1=x+l*st*sp;
      double y1=y-l*st*cp;
      double z1=z+l*ct;  
      double x2=x-l*st*sp;
      double y2=y+l*st*cp;
      double z2=z-l*ct;
      int a = x1*(y2-pointInMu2e.y()) + x2*(pointInMu2e.y()-y1) + pointInMu2e.x()*(y1-y2);
      if(a == 0)
        {
        TEveMu2eCustomHelix *line_twoDstrawXY = new TEveMu2eCustomHelix();
        line_twoDstrawXY->SetLineWidth(1);
        line_twoDstrawXY->SetPoint(0,pointmmTocm(x1),pointmmTocm(y1),pointmmTocm(z1));
        line_twoDstrawXY->SetNextPoint(pointmmTocm(x2),pointmmTocm(y2),pointmmTocm(z2));
        line_twoDstrawXY->SetLineColor(kRed);
        HitList2DXY->AddElement(line_twoDstrawXY);
        
        TEveMu2eCustomHelix *line_twoDstrawXZ = new TEveMu2eCustomHelix();
        line_twoDstrawXZ->SetLineWidth(1);
        line_twoDstrawXZ->SetPoint(0,pointmmTocm(x1),pointmmTocm(y1)+ 1000,pointmmTocm(z1));
        line_twoDstrawXZ->SetNextPoint(pointmmTocm(x2),pointmmTocm(y2)+ 1000,pointmmTocm(z2));
        line_twoDstrawXZ->SetLineColor(kRed);
        HitList2DXZ->AddElement(line_twoDstrawXZ);
        }
      }
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    //this->SetMarkerSize(mSize_);
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
      HitList2DXY->AddElement(error);
    }
    HitList2DXY->AddElement(this);
    HitList2DXZ->AddElement(this);
    }
  }
