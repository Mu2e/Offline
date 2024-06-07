#ifndef TEveMu2eMainWindow_h
#define TEveMu2eMainWindow_h


// ... libRIO
#include <TFile.h>
#include <TObject.h>
#include <TSystem.h>
#include <TText.h>
#include <TCanvas.h>
// ... libGui
#include <TGIcon.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGString.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include <TGLViewer.h>
#include <TGLEmbeddedViewer.h>
#include <TGMsgBox.h>
#include <TGSplitFrame.h>
// ... libRGL
#include <TGLViewer.h>
#include <TVirtualX.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TQObject.h>
// ... libEve
#include <TEvePad.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveParamList.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveStraightLineSet.h>
#include <TEveText.h>
//libGeom
#include <TGeoManager.h>
//Mu2e:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"
//...TEveMu2e

#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCluster.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
#include "Offline/TEveEventDisplay/src/dict_classes/Geom_Interface.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCRV.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eDataInterface.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCInterface.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eProjectionInterface.h"
class TBox;
class TGTextEntry;
class TPad;
class TGCanvas;
class TRootEmbeddedCanvas;
class TGTextButton;
class TText;
class TGSplitFrame;
class TGLOverlayButton;
namespace mu2e{

  struct DrawOptions{
    // data options
    bool addCRVInfo = false;
    bool addCosmicTracks = false;
    bool addTracks = false;
    bool addClusters = false;
    bool addComboHits = false;
    bool addTrkHits = false;
    bool addTimeClusters = false;
    bool addCryHits = false;
    bool addMCTraj = false;
    bool addKKTracks = false;
    DrawOptions(){};
    DrawOptions(bool crv, bool cosmictracks, bool tracks, bool clusters, bool combohits, bool trkhits, bool timeclusters, bool cryhits, bool mctraj, bool kktracks)
    : addCRVInfo(crv), addCosmicTracks(cosmictracks), addTracks(tracks), addClusters(clusters), addComboHits(combohits), addTrkHits(trkhits), addTimeClusters(timeclusters), addCryHits(cryhits), addMCTraj(mctraj), addKKTracks(kktracks) {};
   };



        class TEveMu2eMainWindow : public TGMainFrame {
    public:

      struct GeomOptions{
        GeomOptions(){}
        GeomOptions( bool abuilding, bool aCRV, bool aDSOnly, bool aInsidePS )
          : showbuilding(abuilding)
          , showCRV(aCRV)
          , showDSOnly(aDSOnly)
          , showInsidePS(aInsidePS){
        }

        bool showbuilding = false;
        bool showCRV      = false;
        bool showDSOnly   = false;
        bool showInsidePS = false;
      };


      #ifndef __CINT__
      TEveMu2eMainWindow();
      TEveMu2eMainWindow(const TEveMu2eMainWindow &);
      TEveMu2eMainWindow& operator=(const TEveMu2eMainWindow &);

      TEveMu2eMainWindow(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet _pset, const DrawOptions drawOpts, const GeomOptions geomOpts);
      virtual ~TEveMu2eMainWindow(){};
      enum ETestComandIdentifiers{HId1, HId2, HId3};

      // For viewers:
      void StartProjectionTabs();
      void CreateMultiViews();
      void CreateCaloProjection();
      void CreateTrackerProjection();
      void CreateCRVProjection();
      void PrepareTrackerProjectionTab(const art::Run& run);
      void PrepareCaloProjectionTab(const art::Run& run);
      void PrepareCRVProjectionTab(const art::Run& run);
      void SetParticleOpts(std::vector<int> particles_) { particles = particles_;}

      //GUI and geom:
      void CreateGUI();
      void SetRunGeometry(const art::Run& run, std::string gdmlname, int _diagLevel);
      void RedrawDataProducts(std::string type);

      // for menu:
      Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);

      // to add event info:
      void  setEvent(const art::Event& event, bool firstLoop, Data_Collections &data, double time, bool accumulate, int& runn, int& eventn, bool& update, bool isMCOnly);
      bool  isClosed() const;
      int   getEventToFind(bool &findEvent) const;

      //List of parameters:
      TGeoManager* geom = new TGeoManager("geom","Geom");
      Geom_Interface *mu2e_geom = new Geom_Interface();
      TEveMu2eDataInterface *pass_data  = new TEveMu2eDataInterface();
      TEveMu2eMCInterface *pass_mc      = new TEveMu2eMCInterface();
      TEveMu2eProjectionInterface *pass_proj = new TEveMu2eProjectionInterface();

      // data options
      DrawOptions DrawOpts;
      GeomOptions GeomOpts;
      std::vector<double> *clusterenergy = 0;
      std::vector<double> *hitenergy = 0;
      std::vector<double> times;

      TEvePad *fPad = nullptr;
      TEvePad   *fPadCRV = nullptr;
      TGSplitFrame *fSplitFrame = nullptr;
      TGSplitFrame *fSplitFrameCRV = nullptr;
      TGSplitFrame *frm = nullptr;
      TGSplitFrame *frmCRV = nullptr;
      TGLEmbeddedViewer *fViewer0 = nullptr;
      TGLEmbeddedViewer *fViewer1 = nullptr;
      TGLEmbeddedViewer *fViewer2 = nullptr;
      TGLEmbeddedViewer *fViewer3 = nullptr;
      TGLEmbeddedViewer *fViewer4 = nullptr;
      TGLEmbeddedViewer *fViewer5 = nullptr;
      Bool_t fIsEmbedded;

      TEveViewer *fViewer[6];
      TEveProjectionManager *fRPhiMgr = nullptr;
      TEveProjectionManager *fRhoZMgr = nullptr;
      TEveProjectionManager *fXYMgr = nullptr;
      TEveProjectionManager *gRPhiMgr = 0;
      TEveProjectionManager *gRhoZMgr = 0;
      TEveProjectionManager *TfXYMgr = 0;
      TEveProjectionManager *TfRZMgr = 0;
      TEveProjectionManager *CfXYMgr = 0;
      TEveProjectionManager *CfRZMgr = 0;
      TEveProjectionManager *CrfXYMgr = 0;
      TEveProjectionManager *CrfRZMgr = 0;
      TGMainFrame* frmMain  = nullptr;
      TEveBrowser* browser = nullptr;
      TEveScene *s = 0;
      TEveScene *proj0 = 0;
      TEveScene *proj1 = 0;
      TEveScene *proj2 = 0;
      TEveScene *proj3 = 0;
      TEveScene *proj4 = 0;
      TEveScene *proj5 = 0;


      bool usereventSelected = false;
      TGTextEntry     *fTeRun,*fTeEvt, *fTTEvt, *fTeh1, *fTeh2, *fTeh3, *cminenergy, *cmaxenergy, *hminenergy, *hmaxenergy, *hmintime, *hmaxtime;
      TGLabel         *fTlRun,*fTlEvt, *fTlTEvt, *fTlHSlid, *celabel, *helabel,*timelabel, *spacer, *spacer1;
      TGButtonGroup       *br;
      TGCheckButton       *clusterscheck, *hitscheck, *trackscheck, *cosmicscheck, *cosmictrkscheck, *mctrajcheck;

      TGTextBuffer *_eventNumber, *_subrunNumber, *_runNumber, *_time,  *_clustminenergy, *_clustmaxenergy, *_hitminenergy, *_hitmaxenergy, *_hitmintime, *_hitmaxtime;

      int eventToFind, runToFind;
      int  _eventToFind = 0;

      bool _isClosed = false;
      bool _findEvent = true;
      bool _firstLoop = true;
      bool _accumulate = false;

      TEveMu2e2DProjection *tracker2Dproj = new TEveMu2e2DProjection();
      TEveMu2e2DProjection *calo2Dproj = new TEveMu2e2DProjection();
      TEveMu2e2DProjection *CRV2Dproj = new TEveMu2e2DProjection();

      TEveMu2eCalorimeter *Mu2eCalo = new TEveMu2eCalorimeter();
      TEveMu2eTracker *Mu2eTracker  = new TEveMu2eTracker();
      TEveMu2eCRV *Mu2eCRV = new TEveMu2eCRV();

      int _event, _subrun, _run;
      Data_Collections _data;
      Data_Collections _emptydata;

      double fclustmin = -1;
      double fclustmax = -1;
      double fhitmin = -1;
      double fhitmax = -1;
      double ftimemin = -1;
      double ftimemax = -1;

      std::vector<int> particles;
      #endif
      ClassDef(TEveMu2eMainWindow,0);

    }; //end class def

}//end namespace mu2e

#endif /*TEveMu2eMainWindow.h*/
