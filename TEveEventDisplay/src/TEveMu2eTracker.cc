#include "TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/TargetFoil.hh"

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
      Double_t dr = rmax - rmin;
      
      //Tracker Planes in XZ 
      //int nplane = trkr->getPlane(0).nplanes();
      unsigned int nplanes = trkr->nPlanes();
      double p = 0.0;
      for(int i =0;i<nplanes;i++) 
        { 
        Double_t planepos[3];
        Double_t zplane{pointmmTocm(2*trkr->g4Tracker()->getplaneEnvelopeParams().zHalfLength())};
        TEveGeoShape *plane = new TEveGeoShape();
        CLHEP::Hep3Vector Pos_plane(0,1000,p-dz+zplane);
      
        planepos [0] = Pos_plane.x();
        planepos [1] = Pos_plane.y();
        planepos [2] = Pos_plane.z();
      
        plane->SetShape(new TGeoBBox("plane",rmax+rmin/2,rmax+rmin/2,zplane,planepos));
        plane->SetMainTransparency(100);
        orthodetXZ->AddElement(plane);
        p = p + 15.568;
        }
      
      //XY:
      TEveGeoShape *tr = new TEveGeoShape();
      tr->SetShape(new TGeoTube(rmin, rmax, dz));
      tr->SetMainTransparency(100);
      orthodetXY->AddElement(tr);
      
      //upper
      Double_t origin_upper[3];
      TEveGeoShape *upper = new TEveGeoShape();
      CLHEP::Hep3Vector Pos_upper(0,1000+rmin+dr/2,0);

      origin_upper [0] = Pos_upper.x();
      origin_upper [1] = Pos_upper.y();
      origin_upper [2] = Pos_upper.z();
        
      upper->SetShape(new TGeoBBox("upper",dr,dr,dz,origin_upper));
      upper->SetMainTransparency(100);
      orthodetXZ->AddElement(upper);
      
      //lower
      Double_t origin_lower[3];
      TEveGeoShape *lower = new TEveGeoShape();
      CLHEP::Hep3Vector Pos_lower(0,1000-rmin-dr/2,0);

      origin_lower [0] = Pos_lower.x();
      origin_lower [1] = Pos_lower.y();
      origin_lower [2] = Pos_lower.z();
        
      lower->SetShape(new TGeoBBox("lower",dr,dr,dz,origin_lower));
      lower->SetMainTransparency(100);
      orthodetXZ->AddElement(lower);
        
      // ... Create tracker using the composite shape defined above
      TGeoMaterial *mat = new TGeoMaterial("Mylar", 12,6,1.4);
      TGeoMedium *My = new TGeoMedium("Mylar",2, mat);
      CLHEP::Hep3Vector trackerCentrMu2e = GetTrackerCenter();
      TGeoShape *gs = new TGeoTube("Straw Tracker",rmin,rmax,dz); 
      TGeoVolume *tracker = new TGeoVolume("straw Tracker ",gs, My);
      tracker->SetVisLeaves(kFALSE);
      tracker->SetInvisible();
      topvol->AddNode(tracker, 1, new TGeoTranslation(-390.4,+1000,1017.1)); 
      
      // Addition of Stopping Target geometry
      GeomHandle<StoppingTarget> target;
      CLHEP::Hep3Vector _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();
      double stoppingtargetlength=target->cylinderLength();
      double stoppingtargetz = target->centerInMu2e().z() - _detSysOrigin.z();
      double startz = stoppingtargetz - stoppingtargetlength*0.5;
      unsigned int n=target->nFoils();
      double j =0.0; 
      for(unsigned int i=0; i<n; i++)
        {
        if(i > 0) j = j+abs(target->foil(i-1).centerInMu2e().z() - target->foil(i).centerInMu2e().z());
        const mu2e::TargetFoil &foil=target->foil(i);
        double halfThickness = foil.halfThickness();
        double r = foil.rOut() - foil.rIn();
        CLHEP::Hep3Vector center = foil.centerInDetectorSystem();
        CLHEP::Hep3Vector foilposition(center.x() ,1000+center.y(),startz/10+j/10); // Stopping Target Location 
        Double_t foilpos[3];
        foilpos [0] = foilposition.x();
        foilpos [1] = foilposition.y();
        foilpos [2] = foilposition.z();
      
        TEveGeoShape *stXZ = new TEveGeoShape();
        stXZ->SetShape(new TGeoBBox("foil",pointmmTocm(r),pointmmTocm(r),pointmmTocm(halfThickness), foilpos));
        stXZ->SetMainTransparency(100);
        orthodetXZ->AddElement(stXZ);
      
        TEveGeoShape *stXY = new TEveGeoShape();
        stXY->SetShape(new TGeoTube(pointmmTocm(foil.rIn()),pointmmTocm(foil.rOut()),pointmmTocm(halfThickness)));
        stXY->SetMainTransparency(100);
        orthodetXY->AddElement(stXY);

        
        }


  }
}
