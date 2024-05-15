#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCRV.h"
#include "art/Framework/Principal/Run.h"

using namespace mu2e;
namespace mu2e{
        TEveMu2eCRV::TEveMu2eCRV(){}

  /*------------Function to draw CRV geometry in 2D:-------------*/
  void TEveMu2eCRV::DrawCRVDetector(art::Run const& run, TGeoVolume* topvol, TEveElementList *orthodetT1, TEveElementList *orthodetT2){
    TGeoMaterial *matSi = new TGeoMaterial("Si", 28.085,14,2.33);
    TGeoMedium *Si = new TGeoMedium("Silicon",2, matSi);

    std::vector<double> halflen;
    CLHEP::Hep3Vector position;
    CosmicRayShield const &CRS = *(GeomHandle<CosmicRayShield>());
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
      topvol->AddNode(crv0, 1, new TGeoTranslation(pointmmTocm(Center[0]),pointmmTocm(Center[1]),pointmmTocm(Center[2])));
    }
  }
}

