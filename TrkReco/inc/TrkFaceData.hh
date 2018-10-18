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
    // float  weightXY;
    // float  weightZPhi;
    HitInfo_t(){
      face          = -1;	 
      panel         = -1;
      panelHitIndex = -1;
      // weightXY      = 0.;
      // weightZPhi    = 0.;
    }
    HitInfo_t(int F, int P, int H, float WXY=0., float WZPhi=0.){
      face          = F;	 
      panel         = P;
      panelHitIndex = H;
      // weightXY      = WXY;
      // weightZPhi    = WZPhi;
    }
  };
  
 
  struct PanelZ_t {
    //    int                              fNHits;      // total number of hits per panel
    //    constexpr static uint16_t        kNMaxPanelHits = 20;//maximum number of hits within a panel

    //std::vector<ComboHit>        fHitData;
    int                              idChBegin;
    int                              idChEnd;
    
    int   nChHits(){ return (idChEnd - idChBegin); }

    //    float                            phi;         // phi angle of the wire

    PanelZ_t    (){
      idChBegin = -1;
      idChEnd   = -1;
      //      fNHits  = 0;
      //      fHitData.reserve(kNMaxPanelHits);
    }
  }; 

  struct FaceZ_t {
    constexpr static uint16_t kNPanels           = 3; // number of panels per plane
    constexpr static uint16_t kNPlanesPerStation = 2;

    // float                          z;           // 
    int                       bestFaceHit;             
    std::array<PanelZ_t,kNPanels>  panelZs;
    
    int                              idChBegin;
    int                              idChEnd;
 
    int   nChHits(){ return (idChEnd - idChBegin); }

    FaceZ_t    (){
      bestFaceHit = -1;
      idChBegin   = -1;
      idChEnd     = -1;
    }
    
    int   evalUniqueHitIndex(int &Face, int& Panel, int& PanelHitIndex){
      //      return Face*FaceZ_t::kNPanels*PanelZ_t::kNMaxPanelHits + Panel*PanelZ_t::kNMaxPanelHits + PanelHitIndex;
      return  PanelHitIndex;
    }

    int   evalUniqueHitIndex(HitInfo_t &HitInfo){
      // return HitInfo.face*FaceZ_t::kNPanels*PanelZ_t::kNMaxPanelHits + HitInfo.panel*PanelZ_t::kNMaxPanelHits + HitInfo.panelHitIndex;
      return HitInfo.panelHitIndex;
    }
  }; 


};
#endif
