#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Principal/Run.h"

using namespace mu2e;
namespace mu2e{

    TEveMu2eCalorimeter::TEveMu2eCalorimeter(){}

     /*------------Function to construct calorimeter (for 2D only):-------------*/
    void TEveMu2eCalorimeter::DrawCaloDetector(art::Run const& run, TGeoVolume* topvol, TEveElementList *orthodet0, TEveElementList *orthodet1){
      TGeoMaterial *mat = new TGeoMaterial("CsI", 28.085,14,2.33);
      TGeoMedium *CsI = new TGeoMedium("CsI",2, mat);

      Calorimeter const &cal = *(GeomHandle<Calorimeter>());
      unsigned int CrystalsperDisk = cal.disk(0).nCrystals();
      for(unsigned int i = 0; i < CrystalsperDisk * 2 ; i++){
        Crystal const &crystal = cal.crystal(i);
        double crystalXLen = pointmmTocm(crystal.size().x());
        double crystalYLen = pointmmTocm(crystal.size().y());
        double crystalZLen = pointmmTocm(crystal.size().z());

        CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(0,crystal.position());
        Double_t origin[3];
        GeomHandle<DetectorSystem> det;
        crystalPos = det->toMu2e(crystalPos);

        hep3vectormmTocm(crystalPos);
        origin [0] = crystalPos.x();
        origin [1] = crystalPos.y();
        origin [2] = crystalPos.z();
        TEveGeoShape *crystalShape   = new TEveGeoShape();
        crystalShape->SetMainTransparency(100);
        crystalShape->SetShape(new TGeoBBox("crystalD0", (crystalXLen/2), (crystalYLen/2), (crystalZLen/2)/10, origin));
        if(i < CrystalsperDisk){
          orthodet0->AddElement(crystalShape);
          TGeoShape *c = new TGeoBBox("crystalD0", (crystalXLen/2), (crystalYLen/2), (crystalZLen/2));
          TGeoVolume *cry= new TGeoVolume("cryD0",c, CsI);
          cry->SetVisLeaves(kFALSE);
          cry->SetInvisible();
          topvol->AddNode(cry, 1, new TGeoTranslation(crystalPos.x(),crystalPos.y(),crystalPos.z()));
        } else {
          crystalShape->SetShape(new TGeoBBox("crystalD1", (crystalXLen/2), (crystalYLen/2), (crystalZLen/2), origin));
          orthodet1->AddElement(crystalShape);
          TGeoShape *c = new TGeoBBox("crystalD1", (crystalXLen/2), (crystalYLen/2), (crystalZLen/2));
          TGeoVolume *cry= new TGeoVolume("cryD1",c, CsI);
          cry->SetVisLeaves(kFALSE);
          cry->SetInvisible();
          topvol->AddNode(cry, 1, new TGeoTranslation(crystalPos.x(),crystalPos.y(),crystalPos.z()));
        }
      }
  }

}
