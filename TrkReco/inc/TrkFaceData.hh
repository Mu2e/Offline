//
#ifndef TrkFaceData_HH
#define TrkFaceData_HH

#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include <array>
#include <ostream>
#include <string>
#include <math.h>

namespace mu2e {

  struct HitInfo_t {
    int    face;	
    int    panel;
    int    panelHitIndex;
    float  weightXY;
    float  weightZPhi;
    HitInfo_t(){
      face          = -1;	 
      panel         = -1;
      panelHitIndex = -1;
      weightXY      = 0.;
      weightZPhi    = 0.;
    }
    HitInfo_t(int F, int P, int H, float WXY=0., float WZPhi=0.){
      face          = F;	 
      panel         = P;
      panelHitIndex = H;
      weightXY      = WXY;
      weightZPhi    = WZPhi;
    }
  };
  
 
  struct PanelZ_t {
    int                              fNHits;      // guess, total number of hits per panel
    //    std::vector<ComboHit>            fHitData;
    constexpr static uint16_t        kNMaxPanelHits = 100;//maximum number of hits within a panel

    //    std::array<ComboHit, kNMaxPanelHits>        fHitData;
    std::vector<ComboHit>        fHitData;

    double                           wx;          // direction cosines of the wires, assumed to be all the same
    double                           wy;      
    double                           phi;         // phi angle of the wire
    //    double                           z;           // 

    PanelZ_t    (){
      fNHits  = 0;
      //      fHitData.reserve(kNMaxPanelHits);
    }

    PanelZ_t    (const PanelZ_t&Copy){
      fNHits = Copy.fNHits;  
      //      fPanel = Copy.fPanel;  
      wx     = Copy.wx;      
      wy     = Copy.wy;      
      phi    = Copy.phi;     
      //      z      = Copy.z;       
      int nhits = Copy.fHitData.size();
      //      fHitData.reserve(kNMaxPanelHits);
      for (int i=0; i<nhits; ++i){
	fHitData[i] = ComboHit(Copy.fHitData[i]);
      }
    }
  }; 

  struct FaceZ_t {
    constexpr static uint16_t kNPanels           = 3; // number of panels per plane
    constexpr static uint16_t kNTotalFaces       = StrawId::_nfaces*StrawId::_nplanes;
    constexpr static uint16_t kNPlanesPerStation = 2;

    double                         z;           // 
    HitInfo_t                      bestFaceHit;             
    std::array<PanelZ_t,kNPanels>  panelZs;

    FaceZ_t    (){
      bestFaceHit = HitInfo_t();
    }

    FaceZ_t    (const FaceZ_t&Copy){
      //      fNHits = Copy.fNHits;  
      z           = Copy.z;  
      bestFaceHit = Copy.bestFaceHit;
      for (int i=0; i<kNPanels; ++i){
	panelZs[i] = Copy.panelZs[i];
      }
    }
    
    
    int   evalUniqueHitIndex(int &Face, int& Panel, int& PanelHitIndex){
      return Face*FaceZ_t::kNPanels*PanelZ_t::kNMaxPanelHits + Panel*PanelZ_t::kNMaxPanelHits + PanelHitIndex;
    }

    int   evalUniqueHitIndex(HitInfo_t &HitInfo){
      return HitInfo.face*FaceZ_t::kNPanels*PanelZ_t::kNMaxPanelHits + HitInfo.panel*PanelZ_t::kNMaxPanelHits + HitInfo.panelHitIndex;
    }
  }; 


};
#endif
