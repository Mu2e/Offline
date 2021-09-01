include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
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
  std::tuple<CLHEP::Hep3Vector, CLHEP::Hep3Vector, std::string, int>
TEveMu2eCRVEvent::DrawSciBar(){
         CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
        std::string strawtitle;
        int colorid = 0;
        double dx=0.0,dy=0.0,dz=0.0;
        std::string shieldside = "CRV_D";
          CLHEP::Hep3Vector _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();
          GeomHandle<CosmicRayShield> CRS;
          std::vector<mu2e::CRSScintillatorShield> const& shields = CRS->getCRSScintillatorShields();
    for(std::vector<mu2e::CRSScintillatorShield>::const_iterator ishield=shields.begin(); ishield!=shields.end(); ++ishield)
    {
      CRSScintillatorShield const& shield = *ishield;
      std::string const& shieldName = shield.getName();
       if(shieldName.compare(4,1, shieldside, 4,1) == 0){
      CRSScintillatorBarDetail const& barDetail = shield.getCRSScintillatorBarDetail();
      dx=(barDetail.getHalfLengths()[0]);
      dy=(barDetail.getHalfLengths()[1]);
      dz=(barDetail.getHalfLengths()[2]);
     // std::cout<<dx<<" "<<dy<<" "<<dz<<std::endl;
      }
    }
          Double_t sibarpos[3];
         // const CrvRecoPulseCollection *crvcoincol;
          const CRSScintillatorBarIndex &crvBarIndexn = fCrvRecoPulse_.GetScintillatorBarIndex();
          const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndexn);
          CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
           // GeomHandle<DetectorSystem> det;
	   
	    CLHEP::Hep3Vector barOffset = crvCounterPos - _detSysOrigin;
            sibarpos[0]=(barOffset.x());
            sibarpos[1]=(barOffset.y());// + 1000.0;
            sibarpos[2]=(barOffset.z());  
            int index = 20;
                  strawtitle =Form("index %i", index);   
            colorid = index;
            sposi.set(pointmmTocm(sibarpos[0]-dx), pointmmTocm(sibarpos[1]-dy), pointmmTocm(sibarpos[2]-dz));
            sposf.set(pointmmTocm(sibarpos[0]+dx), pointmmTocm(sibarpos[1]+dy), pointmmTocm(sibarpos[2]+dz));
      
          // }
       
          return {sposi, sposf, strawtitle, colorid};
  }
  /*------------Function to 3D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit3D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList3D)
  {
          auto [sposi, sposf, title, colorid] = DrawSciBar();
    if(sposi.x()!=0){
      //GeomHandle<DetectorSystem> det;
      //CLHEP::Hep3Vector sposin = det->toMu2e(sposi);
      //CLHEP::Hep3Vector sposfn = det->toMu2e(sposf);
      TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
      line->SetLineWidth(1);
      line->SetPoint(0,sposi.x(),sposi.y(),sposi.z());
      line->SetNextPoint(sposf.x(),sposf.y(),sposf.z());
      line->SetLineColor(colorid);
      line->SetTitle(Form(title.c_str()));
      CrvList3D->AddElement(line);
    }
            
    this->SetTitle((DataTitle(pstr, n)).c_str());
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList3D->AddElement(this);
  }

	 /*------------Function to 2D draw hits:-------------*/
  void TEveMu2eCRVEvent::DrawHit2D(const std::string &pstr, Int_t n, CLHEP::Hep3Vector pointInMu2e, TEveElementList *CrvList2DXY, TEveElementList
*CrvList2DXZ)
  {
    auto [sposi, sposf, title, colorid] = DrawSciBar();
    if(sposi.x()!=0){
      TEveMu2eCustomHelix *line_twoDstrawXY = new TEveMu2eCustomHelix();
      line_twoDstrawXY->SetLineWidth(1);
      line_twoDstrawXY->SetPoint(0,(sposi.x()),(sposi.y()),(sposi.z()));
      line_twoDstrawXY->SetNextPoint((sposf.x()),(sposf.y()),(sposf.z()));
      line_twoDstrawXY->SetLineColor(colorid);
      line_twoDstrawXY->SetTitle(Form(title.c_str()));
      CrvList2DXY->AddElement(line_twoDstrawXY);
      
      TEveMu2eCustomHelix *line_twoDstrawXZ = new TEveMu2eCustomHelix();
      line_twoDstrawXZ->SetLineWidth(1);
      line_twoDstrawXZ->SetPoint(0,(sposi.x()),(sposi.y())+1000.0,(sposi.z()));
      line_twoDstrawXZ->SetNextPoint((sposf.x()),(sposf.y())+1000.0,(sposf.z()));
      line_twoDstrawXZ->SetLineColor(colorid);
      line_twoDstrawXZ->SetTitle(Form(title.c_str()));
      CrvList2DXZ->AddElement(line_twoDstrawXZ);
    }
     
     this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z());
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    CrvList2DXY->AddElement(this);
    CrvList2DXZ->AddElement(this);
	  
  }
}

