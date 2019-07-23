
#ifndef ValSimParticleTimeMap_HH_
#define ValSimParticleTimeMap_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValSimParticleTimeMap {

  public:
    ValSimParticleTimeMap(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const SimParticleTimeMap & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _ht;
    TH1D* _ht2;
  };
}


#endif
