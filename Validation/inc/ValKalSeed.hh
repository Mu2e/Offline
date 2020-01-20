#ifndef ValKalSeed_HH_
#define ValKalSeed_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValKalSeed {

  public:
    ValKalSeed(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const KalSeedCollection & coll, art::Event const& event);
    double mcTrkP(art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hNStraw;
    TH1D* _hNSeg;
    TH1D* _hStatus;
    TH1D* _hflt0;
    TH1D* _ht0;
    TH1D* _hchi2;
    TH1D* _hhasCal;
    TH1D* _hfitCon;
    TH1D* _hfitConC;
    TH1D* _hfitConT;
    TH1D* _hp;
    TH1D* _hpC;
    TH1D* _hpT;
    TH1D* _hpce;
    TH1D* _hpcep;
    TH1D* _hpe;
    TH1D* _hD0;
    TH1D* _hPhi0;
    TH1D* _hOmega;
    TH1D* _hZ0;
    TH1D* _hTan;
    TH1D* _hCuts;
    TH1D* _hPRes;
    TH1D* _hCCdisk;
    TH1D* _hCCEoverP;
    TH1D* _hCCDt0;
    TH1D* _hCCDt;
    TH1D* _hCCDOCA;
    TH1D* _hCChlen;
    TH1D* _hCCtlen;
    TH1D* _hHDrift;
    TH1D* _hHDOCA;
    TH1D* _hHEDep;
    TH1D* _hHPanel;
    TH1D* _hSRadLen;
    TH1D* _hSRadLenSum;
  };
}
#endif
