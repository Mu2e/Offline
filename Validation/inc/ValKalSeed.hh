#ifndef ValKalSeed_HH_
#define ValKalSeed_HH_

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValKalSeed {
 public:
  ValKalSeed(std::string name);
  int declare(const art::TFileDirectory& tfs);
  int fill(const KalSeedCollection& coll, art::Event const& event);
  double mcTrkP(art::Event const& event,VirtualDetectorId const& vdid,double& p_pri);
  std::string& name() { return _name; }

 private:
  std::string _name;
  std::map<VirtualDetectorId,SurfaceId> _vdmap;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hNStraw;
  TH1D* _hNSeg;
  TH1D* _hNInter;
  TH1D* _hTraj;
  TH1D* _hStatus;
  TH1D* _ht0;
  TH1D* _ht0e;
  TH1D* _ht02;
  TH1D* _hchi2;
  TH1D* _hhasCal;
  TH1D* _hactiveCal;
  TH1D* _hfitCon;
  TH1D* _hfitConC;
  TH1D* _hfitConT;
  TH1D* _hp;
  TH1D* _hp2;
  TH1D* _hpC;
  TH1D* _hpT;
  TH1D* _hpce;
  TH1D* _hpcep;
  TH1D* _hsignedp;
  TH1D* _hsignedp2;
  TH1D* _hpe;
  TH1D* _hRho;
  TH1D* _hPhi;
  TH1D* _hCost;
  TH1D* _hCuts;
  TH1D* _hPRes;
  TH1D* _hPResA;
  TH1D* _hCCdisk;
  TH1D* _hCCEoverP;
  TH1D* _hCCDt;
  TH1D* _hCCDOCA;
  TH1D* _hCCcdepth;
  TH1D* _hCCtz;
  TH1D* _hHDrift;
  TH1D* _hHDOCA;
  TH1D* _hHEDep;
  TH1D* _hHPanel;
  TH1D* _hSRadLen;
  TH1D* _hSRadLenSum;
};
}  // namespace mu2e
#endif
