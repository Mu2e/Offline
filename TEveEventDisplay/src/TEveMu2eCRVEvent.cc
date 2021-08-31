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
	  
	  CosmicRayShield const &CRS = *(GeomHandle<CosmicRayShield>());  
	  CRSScintillatorBarDetail const& barDetail = shield.getCRSScintillatorBarDetail();
      double dx=barDetail.getHalfLengths()[0];
      double dy=barDetail.getHalfLengths()[1];
      double dz=barDetail.getHalfLengths()[2];
	  
	 std::vector<mu2e::CRSScintillatorShield> const& shields = CRS.getCRSScintillatorShields();
    for(std::vector<mu2e::CRSScintillatorShield>::const_iterator ishield=shields.begin(); ishield!=shields.end(); ++ishield)
    {
      CRSScintillatorShield const& shield = *ishield;
      std::string const& shieldName = shield.getName();
      int nModules = shield.nModules();
      for (int im = 0; im < nModules; ++im)
      {CRSScintillatorModule const & module = shield.getModule(im);
        int nLayers = module.nLayers();
        for (int il = 0; il < nLayers; ++il)
        {CRSScintillatorLayer const & layer = module.getLayer(il);
         int nBars = layer.nBars();
         for (int ib = 0; ib < nBars; ++ib)
         {Double_t sibarpos[3];
           CRSScintillatorBar const & bar = layer.getBar(ib);
           int index = bar.index().asInt();
	   for(unsigned int i=0; i <crvcoincol->size(); i++)
           {
            const CrvRecoPulse &crvRecoPulse = crvcoincol->at(i);
           const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse.GetScintillatorBarIndex();
           if(index== crvBarIndex)
           {std::cout<<"index = "<<index<<std::endl;
            CLHEP::Hep3Vector barOffset = bar.getPosition() - _detSysOrigin;
            double sibarpos[0]=barOffset.x();
            double sibarpos[1]=barOffset.y()+1000.0;
            double sibarpos[2]=barOffset.z();
		  strawtitle =Form("CRV Scintillator %s  module %i  layer %i  bar %i  (index %i)",shieldName.c_str(),im,il,ib, index);
            colorid = index;
            sposi.set(sibarpos[0]-dx, sibarpos[1]-dy, sibarpos[2]-dz);
            sposf.set(sibarpos[0]+dx, sibarpos[1]+dy, sibarpos[2]+dz);
           }
	   
           }
	 }
       }  
      }
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
    auto [sposi, sposf, title, colorid] = DrawSciBar();
    if(sposi.x()!=0){
      TEveMu2eCustomHelix *line_twoDstrawXY = new TEveMu2eCustomHelix();
      line_twoDstrawXY->SetLineWidth(1);
      line_twoDstrawXY->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y()),pointmmTocm(sposi.z()));
      line_twoDstrawXY->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y()),pointmmTocm(sposf.z()));
      line_twoDstrawXY->SetLineColor(colorid);
      line_twoDstrawXY->SetTitle(Form(title.c_str()));
      HitList->AddElement(line_twoDstrawXY);
      
      TEveMu2eCustomHelix *line_twoDstrawXZ = new TEveMu2eCustomHelix();
      line_twoDstrawXZ->SetLineWidth(1);
      line_twoDstrawXZ->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y())+ 1000,pointmmTocm(sposi.z()));
      line_twoDstrawXZ->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y())+ 1000,pointmmTocm(sposf.z()));
      line_twoDstrawXZ->SetLineColor(colorid);
      line_twoDstrawXZ->SetTitle(Form(title.c_str()));
      HitList->AddElement(line_twoDstrawXZ);
    }
	  
	  this->SetTitle((DataTitle(pstr, n)).c_str());
    hep3vectorTocm(pointInMu2e);
    this->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
    this->SetMarkerColor(mColor_);
    this->SetMarkerSize(mSize_);
    this->SetPickable(kTRUE);
    HitList->AddElement(this);
  }
}

