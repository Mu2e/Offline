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

  /*------------Function to display straws which are hit-------*/
  std::tuple<std::array<double, 3>, CLHEP::Hep3Vector, CLHEP::Hep3Vector>TEveMu2eCRVEvent::DrawSciBar(){
  std::array<double, 3> sibarpos;
  std::array<double, 3> sibardet;
  CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
  GeomHandle<CosmicRayShield> CRS;
  const CRSScintillatorBarIndex &crvBarIndexn = fCrvRecoPulse_.GetScintillatorBarIndex();
  const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndexn);
  CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
  const CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
  sibardet[0]=(barDetail.getHalfLengths()[0]);
  sibardet[1]=(barDetail.getHalfLengths()[1]);
  sibardet[2]=(barDetail.getHalfLengths()[2]);
  sibarpos[0]=(crvCounterPos.x());
  sibarpos[1]=(crvCounterPos.y());
  sibarpos[2]=(crvCounterPos.z());
    
  sposi.set((sibarpos[0]-sibardet[0]),(sibarpos[1]-sibardet[1]),(sibarpos[2]-sibardet[2]));
  sposf.set((sibarpos[0]+sibardet[0]),(sibarpos[1]+sibardet[1]),(sibarpos[2]+sibardet[2]));
  return {sibardet, sposi, sposf};
  }

  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList3D)
  {
    auto [sibardet, sposi, sposf] = DrawSciBar();     
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector sposin = det->toMu2e(sposi);
    CLHEP::Hep3Vector sposfn = det->toMu2e(sposf);
    TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
    line->SetLineWidth(1);
    line->SetPoint(0,sposin.x(),sposin.y(),sposin.z());
    line->SetNextPoint(sposfn.x(),sposfn.y(),sposfn.z());
    CrvList3D->AddElement(line);

    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList3D->AddElement(this);
  }

         /*------------Function to 2D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList2DXY, TEveElementList *CrvList2DYZ)
  {
    auto [sibardet, sposi, sposf] = DrawSciBar();  
    Double_t sibarposition[3];
    sibarposition[0] = pointmmTocm(sposi.x()+sposf.x())/2;
    sibarposition[1] = pointmmTocm(sposi.y()+sposf.y())/2;//+1000;
    sibarposition[2] = pointmmTocm(sposi.z()+sposf.z())/2;

    TEveGeoShape *sibar = new TEveGeoShape();
    sibar->SetShape(new
    TGeoBBox("sibar",pointmmTocm(sibardet[0]),pointmmTocm(sibardet[1]),pointmmTocm(sibardet[2]/8), sibarposition));
    sibar->SetMainTransparency(100);
    CrvList2DXY->AddElement(sibar);
    CrvList2DYZ->AddElement(sibar);

    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList2DXY->AddElement(this);
    CrvList2DYZ->AddElement(this);

  }
}



