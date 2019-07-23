
#ifndef ValCrvDigi_HH_
#define ValCrvDigi_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCrvDigi {

  public:
    ValCrvDigi(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CrvDigiCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    TH1D* _hI;
    TH1D* _hIS;
    TH1D* _ht;
    TH1D* _hA;
  };
}


#endif
