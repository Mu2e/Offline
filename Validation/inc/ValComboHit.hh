
#ifndef ValComboHit_HH_
#define ValComboHit_HH_

#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "Validation/inc/ValId.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValComboHit {

  public:
    ValComboHit(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const ComboHitCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    TH1D* _hsid;
    TH1D* _hNcmb;
    TH1D* _hNstr;
    TH1D* _hx;
    TH1D* _hy;
    TH1D* _hz;
    TH1D* _ht;
    TH1D* _hE;
    TH1D* _hqual;
    TH1D* _hwres;
    TH1D* _htres;

  };
}


#endif
