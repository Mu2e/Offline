#ifndef TEveMu2eMainWindow_h
#define TEveMu2eMainWindow_h

#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TText.h>
#include <TGScrollBar.h>
#include <TGSlider.h>
//#include <TGDoubleSlider.h>
//#include <TSlider.h>
//#include <TSliderBox.h>
#include <TCanvas.h>
#include <TQObject.h>
//libGeom
#include <TGeoManager.h>
#include <TBox.h>
#include <TGeoBBox.h>
//TEve
#include <TEveTrack.h>
#include <TEveManager.h>
//fcl:
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
//Mu2e:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
//TEveMu2e
#include "TEveEventDisplay/src/dict_classes/Geom_Interface.h"
#include "TEveEventDisplay/src/dict_classes/Collection_Filler.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eCRV.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eDataInterface.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCInterface.h"
class TBox;
class TGTextEntry;
class TPad;
class TGCanvas;
class TRootEmbeddedCanvas;
class TGTextButton;
class TText;
namespace mu2e{
	class TEveMu2eMainWindow : public TGMainFrame {
    public:
      
      #ifndef __CINT__
      TEveMu2eMainWindow();
      TEveMu2eMainWindow(const TEveMu2eMainWindow &);
      TEveMu2eMainWindow& operator=(const TEveMu2eMainWindow &);
      TEveMu2eMainWindow(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet _pset);
      virtual ~TEveMu2eMainWindow(){};
      enum ETestComandIdentifiers{HId1, HId2, HId3};
      void StartTrackerProjectionTab();
      void PrepareTrackerProjectionTab(const art::Run& run);
      void StartCaloProjectionTab();
      void PrepareCaloProjectionTab(const art::Run& run);
      void StartCRVProjectionTab();
      void PrepareCRVProjectionTab(const art::Run& run);

      void SetRunGeometry(const art::Run& run, int _diagLevel, bool _showBuilding, bool _showDSOnly, bool _showCRV);
      void RedrawDataProducts(std::string type);
      void RedrawGeometry();
      Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);
      void  setEvent(const art::Event& event, bool firstLoop, Data_Collections &data, double time, bool show2D);
      void  fillEvent(bool firstLoop=false);
      bool  isClosed() const;
      int   getEventToFind(bool &findEvent) const;
      double texttime = -1;
      vector<double> *clusterenergy = 0;
      vector<double> *hitenergy = 0;
      #endif

      TGeoManager* geom = new TGeoManager("geom","Geom");
      Geom_Interface *mu2e_geom	=new Geom_Interface(); 
      TEveMu2eDataInterface *pass_data	=new TEveMu2eDataInterface(); 
      TEveMu2eMCInterface *pass_mc	=new TEveMu2eMCInterface(); 
      int eventToFind, runToFind;

      TGTextEntry     *fTeRun,*fTeEvt, *fTTEvt, *fTeh1, *fTeh2, *fTeh3, *cminenergy, *cmaxenergy, *hminenergy, *hmaxenergy;    
      TGHSlider       *fTHSlid;
      TGLabel         *fTlRun,*fTlEvt, *fTlTEvt, *fTlHSlid, *celabel, *helabel, *spacer, *spacer1;
      TGButtonGroup	*br;
      TGCheckButton	*clusterscheck, *hitscheck, *trackscheck, *cosmicscheck;
      Double_t        hitMarkerSize_;
      Double_t        trkMaxR_;
      Double_t        trkMaxZ_;
      Double_t        trkMaxStepSize_;
      Double_t        camRotateCenterH_;
      Double_t        camRotateCenterV_;
      Double_t        camDollyDelta_;
      Int_t	      HSId1;
      TGTextBuffer *_eventNumber, *_subrunNumber, *_runNumber, *_time, *fTbh1, *fTbh2, *fTbh3, *_clustminenergy, *_clustmaxenergy, *_hitminenergy, *_hitmaxenergy;
      int  _eventToFind = 0; ///TODO - this or one above>?

      bool _isClosed = false;
      bool _findEvent = true;
      bool _firstLoop = true;
      bool _show2D = true;
      TEveMu2e2DProjection *tracker2Dproj = new TEveMu2e2DProjection();
      TEveMu2e2DProjection *calo2Dproj = new TEveMu2e2DProjection();
      TEveMu2e2DProjection *CRV2Dproj = new TEveMu2e2DProjection();
      TEveMu2eCalorimeter *Mu2eCalo = new TEveMu2eCalorimeter();
      TEveMu2eTracker *Mu2eTracker  = new TEveMu2eTracker();
      TEveMu2eCRV *Mu2eCRV = new TEveMu2eCRV();

      TText  *_eventNumberText, *_subrunNumberText, *_runNumberText, *_timeText, *_cminenergy;
      int _event, _subrun, _run;
      Data_Collections _data;
      Data_Collections _emptydata;
      std::vector<double> times;
      double fclustmin = -1;
      double fclustmax = -1;
      double fhitmin = -1;
      double fhitmax = -1;


     ClassDef(TEveMu2eMainWindow,0);

    }; //end class def

}//end namespace mu2e

#endif /*TEveMu2eMainWindow.h*/
