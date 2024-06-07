#ifndef Geom_Interface_h
#define Geom_Interface_h

#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGMsgBox.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>

//Geom:
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/GeometryService/inc/WorldG4Maker.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"
#include "Offline/GeometryService/inc/Mu2eHallMaker.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"


namespace mu2e{
        class Geom_Interface {

    public:
      #ifndef __CINT__
      explicit Geom_Interface();
      virtual ~Geom_Interface(){};
      TGeoManager *_geom;
      TGeoManager* getGeom(TString filename) {
        TGeoManager *geom;
        geom = geom->TGeoManager::Import(filename);
        return geom;
      }
      void CreateGeomManager();
      void RemoveComponents();
      void toForeground();
      void InsidePS( TGeoNode * node, bool inPSVac );
      void InsideDS( TGeoNode * node, bool inDSVac );
      void hideTop(TGeoNode* node, int _diagLevel);
      void hideNodesByName(TGeoNode* node, const std::string& str, bool onOff, int _diagLevel) ;
      void showNodesByName(TGeoNode* node, const std::string& str, bool onOff);
      void InsideCRV( TGeoNode * node, bool inCRVVac);
      void hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
      void SolenoidsOnly(TGeoNode* node);
      void TrackerVolumeHeirarchy( TGeoNode * node, std::vector<CLHEP::Hep3Vector> &TransformList );
      CLHEP::Hep3Vector TrackerG4Origin;
      CLHEP::Hep3Vector CaloG4Origin;
      CLHEP::Hep3Vector TrackMu2eOrigin;
      CLHEP::Hep3Vector CaloMu2eOrigin;
      #endif
    ClassDef(Geom_Interface,0);

        }; //end class def

}//end namespace mu2e

#endif /*Geom_Interface.h*/
