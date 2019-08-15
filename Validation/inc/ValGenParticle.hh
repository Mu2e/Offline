
#ifndef ValGenParticle_HH_
#define ValGenParticle_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Validation/inc/ValId.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValGenParticle {

  public:
    ValGenParticle(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const GenParticleCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    ValId _id;
    TH1D* _hp;
    TH1D* _hlogp;
    TH1D* _hx;
    TH1D* _hxt;
    TH1D* _hy;
    TH1D* _hyt;
    TH1D* _hz;
    TH1D* _hzt;

  };
}


#endif
