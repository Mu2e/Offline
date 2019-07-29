
#ifndef ValBkgQual_HH_
#define ValBkgQual_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/BkgQual.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValBkgQual {

  public:
    ValBkgQual(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const BkgQualCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hmva;
    TH1D* _hstat;
  };
}


#endif
