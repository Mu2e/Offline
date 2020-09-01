#include "TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include <TBox.h>
#include <TGeoBBox.h>
using namespace mu2e;
namespace mu2e{

    TEveMu2eTracker::TEveMu2eTracker(){};

    void TEveMu2eTracker::DrawTrackerDetector(art::Run const& run, TGeoVolume* topvol, TEveElementList *orthodet){
      GeomHandle<Tracker> trkr;

      TubsParams envelope(trkr->getInnerTrackerEnvelopeParams());

      TGeoMaterial *mat = new TGeoMaterial("Mylar", 12,6,1.4);
      TGeoMedium *My = new TGeoMedium("Mylar",2, mat);

      Double_t dz{pointmmTocm(envelope.zHalfLength())};
      Double_t rmin{pointmmTocm(envelope.innerRadius())};
      Double_t rmax{pointmmTocm(envelope.outerRadius())};
      TEveGeoShape *tr = new TEveGeoShape();
      tr->SetShape(new TGeoTube(rmin, rmax, dz));
      tr->SetMainTransparency(100);
      orthodet->AddElement(tr);
        
      // ... Create tracker using the composite shape defined above
      CLHEP::Hep3Vector trackerCenterGDML = GetGDMLTrackerCenter();
      CLHEP::Hep3Vector trackerCentrMu2e = GetTrackerCenter();
      TGeoShape *gs = new TGeoTube("Straw Tracker",rmin,rmax,dz); 
      TGeoVolume *tracker = new TGeoVolume("straw Tracker ",gs, My);
      tracker->SetVisLeaves(kFALSE);
      tracker->SetInvisible();
      topvol->AddNode(tracker, 1, new TGeoTranslation(-390.4,0,1017.1)); 

  }
}
