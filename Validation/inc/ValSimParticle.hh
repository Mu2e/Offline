
#ifndef ValSimParticle_HH_
#define ValSimParticle_HH_

#include "art/Framework/Principal/Event.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Validation/inc/ValId.hh"
#include "Validation/inc/ValPosition.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValSimParticle {

  public:
    ValSimParticle(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const SimParticleCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    ValId _id;
    TH1D* _hp;
    ValPosition _pos;
    TH1D* _hscode;
    TH1D* _hecode;
    ValId _idh;
    TH1D* _hscodeh;
    TH1D* _hecodeh;
  };
}


#endif
