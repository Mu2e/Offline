
#ifndef ValCaloCrystalHit_HH_
#define ValCaloCrystalHit_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCaloCrystalHit {

  public:
    ValCaloCrystalHit(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CaloCrystalHitCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hI;
    TH1D* _ht;
    TH1D* _hE;
  };
}


#endif
