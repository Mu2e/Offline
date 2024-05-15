#include "Offline/GeometryService/inc/GeomHandle.hh"
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
  std::tuple<CLHEP::Hep3Vector, CLHEP::Hep3Vector>TEveMu2eCRVEvent::DrawSciBar(){
  GeomHandle<CosmicRayShield> CRS;
  const CRSScintillatorBarIndex &crvBarIndexn = fCrvRecoPulse_.GetScintillatorBarIndex();
  const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndexn);
  CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
  hep3vectormmTocm(crvCounterPos);
  const CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
  CLHEP::Hep3Vector sibardetails(pointmmTocm(barDetail.getHalfLengths()[0]),pointmmTocm(barDetail.getHalfLengths()[1]),pointmmTocm(barDetail.getHalfLengths()[2]));
  return {sibardetails, crvCounterPos};
  }

  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList3D)
  {
    auto [sibardetails, crvCounterPos] = DrawSciBar();
    Double_t sibarposition[3];
    sibarposition[0] = (crvCounterPos.x());
    sibarposition[1] = (crvCounterPos.y());
    sibarposition[2] = (crvCounterPos.z());
    TEveGeoShape *sibar = new TEveGeoShape();
    sibar->SetShape(new TGeoBBox("sibar",sibardetails.x(),sibardetails.y(),sibardetails.z(), sibarposition));
    CrvList3D->AddElement(sibar);

    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList3D->AddElement(this);
  }

         /*------------Function to 2D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit2DXY(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList2DXY)
  {
    auto [sibardetails, crvCounterPos] = DrawSciBar();
    Double_t sibarposition[3];
    sibarposition[0] = pointmmTocm(crvCounterPos.x());
    sibarposition[1] = pointmmTocm(crvCounterPos.y());
    sibarposition[2] = pointmmTocm(crvCounterPos.z());

    TEveGeoShape *sibar = new TEveGeoShape();
    sibar->SetShape(new TGeoBBox("sibar",pointmmTocm(sibardetails.x()),pointmmTocm(sibardetails.y()),pointmmTocm(sibardetails.z()), sibarposition));
    sibar->SetMainTransparency(100);
    CrvList2DXY->AddElement(sibar);

    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectormmTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList2DXY->AddElement(this);
  }

  void TEveMu2eCRVEvent::DrawHit2DYZ(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList2DYZ)
  {
    auto [sibardetails, crvCounterPos] = DrawSciBar();
    Double_t sibarposition[3];
    sibarposition[0] = 0.0;
    sibarposition[1] = pointmmTocm(crvCounterPos.y());
    sibarposition[2] = pointmmTocm(crvCounterPos.z());

    TEveGeoShape *sibar = new TEveGeoShape();
    sibar->SetShape(new TGeoBBox("sibar",0.0,pointmmTocm(sibardetails.y()),pointmmTocm(sibardetails.z()), sibarposition));
    sibar->SetMainTransparency(100);
    CrvList2DYZ->AddElement(sibar);

    this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectormmTocm(pointInMu2e);
    this->SetNextPoint(0.0, pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList2DYZ->AddElement(this);

  }
}
