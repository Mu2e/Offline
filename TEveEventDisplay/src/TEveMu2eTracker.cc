#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
Int_t transpOpt = 100;
using namespace mu2e;
namespace mu2e{

    TEveMu2eTracker::TEveMu2eTracker(){}

     /*------------Function to construct Tracker (for 2D only):-------------*/
    void TEveMu2eTracker::DrawTrackerDetector(TGeoVolume* topvol, TEveElementList *orthodetXZ, TEveElementList *orthodetXY){
      GeomHandle<Tracker> trkr;

      TubsParams envelope(trkr->g4Tracker()->getInnerTrackerEnvelopeParams());
      TubsParams planeEnvelope(trkr->g4Tracker()->getPlaneEnvelopeParams());

      Double_t dz{pointmmTocm(envelope.zHalfLength())};
      Double_t zpanel{pointmmTocm(planeEnvelope.zHalfLength())};
      Double_t rmin{pointmmTocm(planeEnvelope.innerRadius())};
      Double_t rmax{pointmmTocm(planeEnvelope.outerRadius())};

      Double_t fullLength = trkr->planes().back().origin().z()- trkr->planes().front().origin().z();
      Double_t planespace = 2*((pointmmTocm(fullLength) - trkr->nPlanes()*zpanel)/(trkr->nPlanes()-2));
      double p = -dz;
      Double_t panelpos[3];
      panelpos [0] = 0.0;
      panelpos [1] = 0.0;
      //Tracker Planes in XZ
      for(size_t i =0;i<trkr->nPlanes()/2;i++)
        {

        TEveGeoShape *panel = new TEveGeoShape();
        panelpos [2] = p;
        panel->SetShape(new TGeoBBox("panel",rmax+rmin/2,rmax+rmin/2,zpanel*2,panelpos));
        panel->SetMainTransparency(transpOpt);
        orthodetXZ->AddElement(panel);
        p = p + planespace + zpanel*2;
        }

      //XY:
      TEveGeoShape *tr = new TEveGeoShape();
      tr->SetShape(new TGeoTube(rmin, rmax, dz));
      tr->SetMainTransparency(transpOpt);
      orthodetXY->AddElement(tr);

      // ... Create tracker using the composite shape defined above
      TGeoMaterial *mat = new TGeoMaterial("Mylar", 12,6,1.4);
      TGeoMedium *My = new TGeoMedium("Mylar",2, mat);
      GeomHandle<DetectorSystem> det;
      CLHEP::Hep3Vector trackerCentrMu2e = det->getOrigin();

      TGeoShape *gs = new TGeoTube("Straw Tracker",rmin,rmax,dz);
      TGeoVolume *tracker = new TGeoVolume("Straw Tracker ",gs, My);
      tracker->SetVisLeaves(kFALSE);
      tracker->SetInvisible();
      topvol->AddNode(tracker, 1);

      //Stopping Target
      GeomHandle<StoppingTarget> target;
      CLHEP::Hep3Vector _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();
      double stoppingtargetlength = target->cylinderLength();
      double stoppingtargetz = target->centerInMu2e().z() - _detSysOrigin.z();
      double startz = stoppingtargetz - stoppingtargetlength*0.5;
      unsigned int nFoils = target->nFoils();
      double j =0.0;
      for(unsigned int i=0; i<nFoils; i++)
        {
        if(i > 0) j = j+abs(target->foil(i-1).centerInMu2e().z() - target->foil(i).centerInMu2e().z());
        const mu2e::TargetFoil &foil=target->foil(i);
        double halfThickness = foil.halfThickness();
        double r = foil.rOut() - foil.rIn();

        CLHEP::Hep3Vector center = foil.centerInDetectorSystem();

        CLHEP::Hep3Vector foilposition(center.x() ,center.y(),pointmmTocm(startz+j)); // Stopping Target Location

        Double_t foilpos[3];
        foilpos [0] = foilposition.x();
        foilpos [1] = foilposition.y();
        foilpos [2] = foilposition.z();

        TEveGeoShape *stXZ = new TEveGeoShape();

        stXZ->SetShape(new TGeoBBox("foil",pointmmTocm(r),pointmmTocm(r),pointmmTocm(halfThickness), foilpos));
        stXZ->SetMainTransparency(transpOpt);
        orthodetXZ->AddElement(stXZ);

        TEveGeoShape *stXY = new TEveGeoShape();
        stXY->SetShape(new TGeoTube(pointmmTocm(foil.rIn()),pointmmTocm(foil.rOut()),pointmmTocm(halfThickness)));

        stXY->SetMainTransparency(transpOpt);
        orthodetXY->AddElement(stXY);
        }

        //Calo disk outline in the Tracker XZ display window
        GeomHandle<DiskCalorimeter> calo;

        double diskInnerRingIn = calo->caloInfo().getDouble("diskInnerRingIn");
        double diskOuterRingOut = calo->caloInfo().getDouble("diskOuterRingOut");
        double diskOuterRailOut = diskOuterRingOut + calo->caloInfo().getDouble("diskOutRingEdgeRLength");
        double FPCarbonDZ = calo->caloInfo().getDouble("FPCarbonZLength")/2.0;
        double FPFoamDZ = calo->caloInfo().getDouble("FPFoamZLength")/2.0;
        double FPCoolPipeRadius = calo->caloInfo().getDouble("FPCoolPipeRadius");
        double pipeRadius = calo->caloInfo().getDouble("pipeRadius");
        double frontPanelHalfThick = (2.0*FPCarbonDZ+2.0*FPFoamDZ-pipeRadius+FPCoolPipeRadius)/2.0;
        double holeDZ = calo->caloInfo().getDouble("BPHoleZLength")/2.0;

        double crystalDXY = calo->caloInfo().getDouble("crystalXYLength")*2.0;
        double crystalDZ = calo->caloInfo().getDouble("crystalZLength")/2.0;
        double wrapperDXY = crystalDXY + calo->caloInfo().getDouble("wrapperThickness");

        Double_t crystalpos[3];
        for(unsigned int idisk=0; idisk<calo->nDisk(); idisk++){

          CLHEP::Hep3Vector diskPos = calo->disk(idisk).geomInfo().origin() + CLHEP::Hep3Vector(0.0, 0.0, (-holeDZ+frontPanelHalfThick)) - _detSysOrigin;
          double diskXZwidth = diskOuterRailOut + diskInnerRingIn;
          crystalpos [0] = pointmmTocm(diskPos.x());
          crystalpos [2] = pointmmTocm(diskPos.z());
          while(diskXZwidth > -(diskOuterRailOut + diskInnerRingIn))
            {
            crystalpos [1] = pointmmTocm(diskPos.y() + diskXZwidth);

            TEveGeoShape *crystalXZ = new TEveGeoShape();
            crystalXZ->SetShape(new TGeoBBox("Crystal",pointmmTocm(wrapperDXY),pointmmTocm(wrapperDXY),pointmmTocm(crystalDZ), crystalpos));
            crystalXZ->SetMainTransparency(transpOpt);
            orthodetXZ->AddElement(crystalXZ);
            diskXZwidth-=wrapperDXY;
            }
        }
   }
}


