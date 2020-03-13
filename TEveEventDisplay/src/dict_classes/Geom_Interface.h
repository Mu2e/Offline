#ifndef Geom_Interface_h
#define Geom_Interface_h

#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
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

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"



class TGeoManager;
class TGeoVolume;
class TGMainFrame;

namespace mu2e{
	class Geom_Interface {
                #ifndef __CINT__
		explicit Geom_Interface();
		
		public:
		  virtual ~Geom_Interface(){};
		private:
		  art::Event  *_event;
		  art::Run    *_run;
	  	  TGeoManager *_geom;
		 
		  void CreateGeomManager();
		  void RemoveComponents();
		  void toForeground();
		  void InsideDS( TGeoNode * node, bool inDSVac );
		  void hideTop(TGeoNode* node);
		  void hideNodesByName(TGeoNode* node, const std::string& str,bool onOff) ;
		  void hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
		  void hideBuilding(TGeoNode* node);
		  #endif
     		ClassDef(Geom_Interface,0);

	}; //end class def

}//end namespace mu2e

#endif /*Geom_Interface.h*/
