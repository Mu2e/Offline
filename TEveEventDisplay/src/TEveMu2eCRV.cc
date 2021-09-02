#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCRV.h"

using namespace mu2e;
namespace mu2e{
        TEveMu2eCRV::TEveMu2eCRV(){};

  /*------------Function to draw CRV geometry in 2D:-------------*/
  void TEveMu2eCRV::DrawCRVDetector(art::Run const& run, TGeoVolume* topvol, TEveElementList *orthodetT1, TEveElementList *orthodetT2){
    TGeoMaterial *matSi = new TGeoMaterial("Si", 28.085,14,2.33);
    TGeoMedium *Si = new TGeoMedium("Silicon",2, matSi);
   // CLHEP::Hep3Vector _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();
    std::vector<double> halflen;
    CLHEP::Hep3Vector position;
    CosmicRayShield const &CRS = *(GeomHandle<CosmicRayShield>());

           std::vector<mu2e::CRSScintillatorShield> const& shields = CRS.getCRSScintillatorShields();
    for(std::vector<mu2e::CRSScintillatorShield>::const_iterator ishield=shields.begin(); ishield!=shields.end(); ++ishield)
    {
      CRSScintillatorShield const& shield = *ishield;
      std::string const& shieldName = shield.getName();
     // std::cout<<"Shield name = "<<shieldName<<std::endl;
      CRSScintillatorBarDetail const& barDetail = shield.getCRSScintillatorBarDetail();
      double dx=barDetail.getHalfLengths()[0];
      double dy=barDetail.getHalfLengths()[1];
      double dz=barDetail.getHalfLengths()[2];
      std::string shieldside = "CRV_D";
     // std::cout<<"sci bar details :"<<dx<<" "<<dy<<" "<<dz<<std::endl;
      int nModules = shield.nModules();
      for (int im = 0; im < nModules; ++im)
      {
        CRSScintillatorModule const & module = shield.getModule(im);
        int nLayers = module.nLayers();
        for (int il = 0; il < nLayers; ++il)
        {
          CRSScintillatorLayer const & layer = module.getLayer(il);
          int nBars = layer.nBars();
          for (int ib = 0; ib < nBars; ++ib)
          { Double_t sibarpos[3];
            CRSScintillatorBar const & bar = layer.getBar(ib);
            int index = bar.index().asInt();
            CLHEP::Hep3Vector barOffset = bar.getPosition();// - _detSysOrigin;
            sibarpos[0]=barOffset.x();
            sibarpos[1]=barOffset.y();// +1000.0;
            sibarpos[2]=barOffset.z();

            //boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
            std::string c=Form("CRV Scintillator %s  module %i  layer %i  bar %i  (index %i)",shieldName.c_str(),im,il,ib, index);
            TEveGeoShape *sibar = new TEveGeoShape();
            sibar->SetShape(new TGeoBBox("sibar",pointmmTocm(dx),pointmmTocm(dy),pointmmTocm(dz), sibarpos));
        sibar->SetMainTransparency(100);
     
        // if(strcmp(shieldName.c_str(), shieldside.c_str()) > 0) {std::cout<<shieldName<<std::endl;}
         if(shieldName.compare(4,1, shieldside, 4,1) == 0){
          orthodetT2->AddElement(sibar);
          }
          }
        }
      }  
    }
          
    const std::string TopSectorNames[] = {"T1", "T2", "T3", "T4"};
    for (unsigned int i=0; i<4; i++){
            Double_t panelpos[3];
      halflen = CRS.getSectorHalfLengths(TopSectorNames[i]);
      position = CRS.getSectorPosition(TopSectorNames[i]);
            
        panelpos [0] = position.x();
        panelpos [1] = position.y();  
        panelpos [2] = position.z();
            
      TEveGeoShape *sectorshape = new TEveGeoShape();
      sectorshape->SetShape(new TGeoBBox("sectorshape",pointmmTocm(2*halflen[0]), pointmmTocm(2*halflen[2]), pointmmTocm(2*halflen[1]),panelpos));
      sectorshape->SetMainTransparency(100);
      TGeoShape *g = new TGeoBBox("CRV Sector",pointmmTocm(2*halflen[0]), pointmmTocm(2*halflen[2]), pointmmTocm(2*halflen[1]));
      TGeoVolume *crv0= new TGeoVolume("CRV Sector",g, Si);
      crv0->SetVisLeaves(kFALSE);
      crv0->SetInvisible();
      std::string filename("Mu2eG4/geom/crv_counters_v07.txt");
      SimpleConfig Config(filename);
      std::vector<double> Center;
      if(i==0)  Config.getVectorDouble("crs.firstCounterT1", Center);
      if(i==1)  Config.getVectorDouble ("crs.firstCounterT2",Center);
      if(i==2)  Config.getVectorDouble("crs.firstCounterT3", Center);
      if(i==3)  Config.getVectorDouble("crs.firstCounterT4", Center) ;
       orthodetT1->AddElement(sectorshape);
       topvol->AddNode(crv0, 1, new TGeoTranslation(pointmmTocm(Center[0]),pointmmTocm(Center[1]),pointmmTocm(Center[2]/10)));
    }
  }
}
        

