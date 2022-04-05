
#ifndef ValEventWindowMarker_HH_
#define ValEventWindowMarker_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValEventWindowMarker {

  public:
    ValEventWindowMarker(std::string name):_name(name){}
    int declare( const art::TFileDirectory& tfs);
    int fill(const EventWindowMarker & obj, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _hst;
    TH1D* _hlen;

  };
}


#endif
