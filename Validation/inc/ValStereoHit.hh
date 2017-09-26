
#ifndef ValStereoHit_HH_
#define ValStereoHit_HH_

#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "Validation/inc/ValId.hh"
#include "Validation/inc/ValPosition.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValStereoHit {

  public:
    ValStereoHit(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const StereoHitCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    ValPosition _pos;
    TH1D* _ht;
    TH1D* _hE;
    TH1D* _hchi2;
    TH1D* _hmva;

  };
}


#endif
