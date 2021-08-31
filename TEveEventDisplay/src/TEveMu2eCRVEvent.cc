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
  std::tuple<CLHEP::Hep3Vector, CLHEP::Hep3Vector, std::string, int> TEveMu2eHit::DrawSciBar(){
         CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
        std::string strawtitle;
        int colorid = 0;
	  CLHEP::Hep3Vector _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();
	  GeomHandle<CosmicRayShield> CRS;
	  Double_t sibarpos[3];
          for(unsigned int i=0; i <crvcoincol->size(); i++)
           {
            const CrvRecoPulse &crvRecoPulse = crvcoincol->at(i);
           const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse.GetScintillatorBarIndex();
		  const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
		  int index = 1;
            CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
            GeomHandle<DetectorSystem> det;
            CLHEP::Hep3Vector barOffset = crvCounterPos - _detSysOrigin;
            sibarpos[0]=barOffset.x();
            sibarpos[1]=barOffset.y()+1000.0;
            sibarpos[2]=barOffset.z();
		  strawtitle =Form("index %i", index);
            colorid = index;
            sposi.set(sibarpos[0]-20, sibarpos[1]-20, sibarpos[2]-20);
            sposf.set(sibarpos[0]+20, sibarpos[1]+20, sibarpos[2]+20);
           
           }
	
	  return {sposi, sposf, strawtitle, colorid};
  }
  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *HitList)
  {
	  auto [sposi, sposf, title, colorid] = DrawSciBar();
    if(sposi.x()!=0){
      GeomHandle<DetectorSystem> det;
      CLHEP::Hep3Vector sposin = det->toMu2e(sposi);
      CLHEP::Hep3Vector sposfn = det->toMu2e(sposf);
      TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
      line->SetLineWidth(1);
      line->SetPoint(0,sposin.x(),sposin.y(),sposin.z());
      line->SetNextPoint(sposfn.x(),sposfn.y(),sposfn.z());
      line->SetLineColor(colorid);
      line->SetTitle(Form(title.c_str()));
      HitList->AddElement(line);
    }
	    
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

