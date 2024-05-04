#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"

using namespace mu2e;
namespace mu2e{

  TEveMu2eHit::TEveMu2eHit(){}

  /*------------Function to build title:-------------*/
  std::string TEveMu2eHit::DataTitle(const std::string &pstr, int n){
        std::string dstr=" hit#" + std::to_string(n) + "\nLayer: ";
        std::string strlab=pstr+dstr;
        return (strlab);
  }

  /*------------Function to display straws which are hit-------*/
  std::tuple<CLHEP::Hep3Vector, CLHEP::Hep3Vector, std::string, int> TEveMu2eHit::DrawStraw(){
        mu2e::GeomHandle<mu2e::Tracker> tracker;
        const auto& allStraws = tracker->getStraws();
        CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
        std::string strawtitle;
        int colorid = 0;
        for (size_t i = 0; i<tracker->nStraws(); i++){
          const mu2e::Straw& s = allStraws[i];
          if(s.id().asUint16()==fComboHit_._sid.asUint16())
          {
          const CLHEP::Hep3Vector& p = s.getMidPoint();
          const CLHEP::Hep3Vector& d = s.getDirection();
          int idStraw =  s.id().getStraw();
          int idPanel =  s.id().getPanel();
          int idPlane =  s.id().getPlane();
          colorid = idPlane + idPanel;
          strawtitle =Form("Straw %i Panel %i  Plane %i",idStraw,idPanel,idPlane);
          //std::cout<<idPanel<<" "<<idPlane<<std::endl;
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
          sposi.set(x1, y1, z1);
          sposf.set(x2, y2, z2);
          }
        }
        return {sposi, sposf, strawtitle, colorid};
  }

  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eHit::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList)
  {
    auto [sposi, sposf, title, colorid] = DrawStraw();
    if(sposi.x()!=0){
      GeomHandle<DetectorSystem> det;
      CLHEP::Hep3Vector sposin = det->toMu2e(sposi);
      CLHEP::Hep3Vector sposfn = det->toMu2e(sposf);
      TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
      line->SetLineWidth(1);
      line->SetPoint(0,pointmmTocm(sposin.x()),pointmmTocm(sposin.y()),pointmmTocm(sposin.z()));
      line->SetNextPoint(pointmmTocm(sposfn.x()),pointmmTocm(sposfn.y()),pointmmTocm(sposfn.z()));
      line->SetLineColor(colorid);
      line->SetTitle(Form(title.c_str()));
      HitList->AddElement(line);
    }

    this->SetTitle((DataTitle(pstr, n)).c_str());
    //hep3vectormmTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    this->SetPickable(kTRUE);
    if(AddErrorBar_){
      TEveLine *error = new TEveLine();
      auto const& p = fComboHit_.pos();
      auto w = fComboHit_.uDir();
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
 void TEveMu2eHit::DrawHit2DXY(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList2DXY)
  {
    auto [sposi, sposf, title, colorid] = DrawStraw();
    if(sposi.x()!=0){
      TEveMu2eCustomHelix *line_twoDstrawXY = new TEveMu2eCustomHelix();
      line_twoDstrawXY->SetLineWidth(1);
      line_twoDstrawXY->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y()),pointmmTocm(sposi.z()));
      line_twoDstrawXY->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y()),pointmmTocm(sposf.z()));
      line_twoDstrawXY->SetLineColor(colorid);
      line_twoDstrawXY->SetTitle(Form(title.c_str()));
      HitList2DXY->AddElement(line_twoDstrawXY);

    }
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectormmTocm(pointInMu2e);
    this->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()), pointmmTocm(pointInMu2e.z()));
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    this->SetPickable(kTRUE);

    if(AddErrorBar_){
      TEveLine *error = new TEveLine();
      auto const& p = fComboHit_.pos();
      auto w = fComboHit_.uDir();
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
  }

 void TEveMu2eHit::DrawHit2DXZ(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, int energylevel, TEveElementList *HitList2DXZ)
  {
    auto [sposi, sposf, title, colorid] = DrawStraw();
    if(sposi.x()!=0){
      TEveMu2eCustomHelix *line_twoDstrawXZ = new TEveMu2eCustomHelix();
      line_twoDstrawXZ->SetLineWidth(1);
      line_twoDstrawXZ->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y()),pointmmTocm(sposi.z()));
      line_twoDstrawXZ->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y()),pointmmTocm(sposf.z()));
      line_twoDstrawXZ->SetLineColor(colorid);
      line_twoDstrawXZ->SetTitle(Form(title.c_str()));
      HitList2DXZ->AddElement(line_twoDstrawXZ);
    }
    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectormmTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    int colors[] = {-7, 3, -6, -1, 9, 0, -4, 10, 1};
    this->SetMarkerColor(kSpring + colors[energylevel]);
    this->SetPickable(kTRUE);
    HitList2DXZ->AddElement(this);
  }
}

