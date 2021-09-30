#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"

Int_t transpOpt = 100;
using namespace mu2e;
namespace mu2e{

    TEveMu2eTracker::TEveMu2eTracker(){};

     /*------------Function to construct Tracker (for 2D only):-------------*/
    void TEveMu2eTracker::DrawTrackerDetector(TGeoVolume* topvol, TEveElementList *orthodetXZ, TEveElementList *orthodetXY){
      GeomHandle<Tracker> trkr;

      TubsParams envelope(trkr->g4Tracker()->getInnerTrackerEnvelopeParams());

      Double_t dz{pointmmTocm(envelope.zHalfLength())};
      Double_t rmin{pointmmTocm(envelope.innerRadius())};
      Double_t rmax{pointmmTocm(envelope.outerRadius())};
        
      //Tracker Stations in XZ
      double p = 0.0;
      unsigned int nStations = trkr->nPlanes()/2;
      Double_t stationpos[3];
      for(unsigned int i =0; i < nStations ;i++){
        
        Double_t zpanel{pointmmTocm(2*trkr->g4Tracker()->getPanelEnvelopeParams().zHalfLength())};
        TEveGeoShape *station= new TEveGeoShape();
        CLHEP::Hep3Vector Pos_station(0,1000,p-dz+zpanel);
        int dp = Pos_station.z() - stationpos [2];
        stationpos [0] = Pos_station.x();
        stationpos [1] = Pos_station.y();
        stationpos [2] = Pos_station.z();

        station->SetShape(new TGeoBBox("station",rmax+rmin/2,rmax+rmin/2,zpanel,stationpos));
        station->SetMainTransparency(transpOpt); 
        orthodetXZ->AddElement(station);
        p = p + abs(dp)/10; 
      }

        
      //XY:
      TEveGeoShape *tr = new TEveGeoShape();
      tr->SetShape(new TGeoTube(rmin, rmax, dz));
      tr->SetMainTransparency(transpOpt); 
      orthodetXY->AddElement(tr);

      // ... Create tracker using the composite shape defined above
      TGeoMaterial *mat = new TGeoMaterial("Mylar", 12,6,1.4);
      TGeoMedium *My = new TGeoMedium("Mylar",2, mat);
      CLHEP::Hep3Vector trackerCentrMu2e = GetTrackerCenter();
      TGeoShape *gs = new TGeoTube("Straw Tracker",rmin,rmax,dz);
      TGeoVolume *tracker = new TGeoVolume("Straw Tracker ",gs, My);
      tracker->SetVisLeaves(kFALSE);
      tracker->SetInvisible();
      topvol->AddNode(tracker, 1, new TGeoTranslation(-390.4,+1000,1017.1));  // FIXME - hardcoded number
     
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
        CLHEP::Hep3Vector foilposition(center.x() ,1000+center.y(),startz/10+j/10); // Stopping Target Location FIXME
      
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
        
        for(unsigned int idisk=0; idisk<calo->nDisk(); idisk++){
          CLHEP::Hep3Vector diskPos = calo->disk(idisk).geomInfo().origin() - _detSysOrigin;
          diskPos += CLHEP::Hep3Vector(0.0, 10000.0, (-holeDZ+frontPanelHalfThick));
          double k = diskOuterRailOut + diskInnerRingIn;
          for(unsigned int ic=0; ic < 30; ic++){ // FIXME - what is 30?
            CLHEP::Hep3Vector crystalposition(diskPos.x(), (diskPos.y() + k), diskPos.z());
            Double_t crystalpos[3];
            crystalpos [0] = pointmmTocm(crystalposition.x());
            crystalpos [1] = pointmmTocm(crystalposition.y());
            crystalpos [2] = pointmmTocm(crystalposition.z());
        
            TEveGeoShape *crystalXZ = new TEveGeoShape();
            crystalXZ->SetShape(new TGeoBBox("Crystal",pointmmTocm(wrapperDXY),pointmmTocm(wrapperDXY),pointmmTocm(crystalDZ), crystalpos));
            crystalXZ->SetMainTransparency(transpOpt);   
            orthodetXZ->AddElement(crystalXZ);
            k = k - wrapperDXY;
            }
          }
  }
}


