
#ifndef ValHelixSeed_HH_
#define ValHelixSeed_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValHelixSeed {

  public:
    ValHelixSeed(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const HelixSeedCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hNCombo;
    TH1D* _hNStrHit;
    TH1D* _hStatus;
    TH1D* _ht0;
    TH1D* _hp;
    TH1D* _hpce;
    TH1D* _hpt;
    TH1D* _hD0;
    TH1D* _hPhi0;
    TH1D* _hLambda;
    TH1D* _hchi2dXY;
    TH1D* _hchi2dZPhi;
  };
}


#endif
