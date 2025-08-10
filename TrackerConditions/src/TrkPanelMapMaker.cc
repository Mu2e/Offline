// clang-format off
#include "Offline/TrackerConditions/inc/TrkPanelMapEntity.hh"
#include "Offline/TrackerConditions/inc/TrkPanelMapMaker.hh"
// #include "cetlib_except/exception.h"
// #include "TMath.h"
// #include <cmath>
// #include <complex>
#include <memory>
#include <iostream>

using namespace std;

namespace mu2e {

//-----------------------------------------------------------------------------  
// all vectors are supposed to have the same length  
  TrkPanelMapEntity::ptr_t TrkPanelMapMaker::fromFcl() {

    auto ptr = std::make_shared<TrkPanelMapEntity>();

    std::vector<int> mnid    = config_.mnid   ();
    std::vector<int> dtcid   = config_.dtcid  ();
    std::vector<int> link    = config_.link   ();
    std::vector<int> station = config_.station();
    std::vector<int> psid    = config_.psid   ();
    std::vector<int> plane   = config_.plane  ();
    std::vector<int> ppid    = config_.ppid   ();
    std::vector<int> panel   = config_.panel  ();
    std::vector<int> zface   = config_.zface  ();
    
    int npanels = mnid.size();
    for (int i=0; i<npanels; i++) {
      TrkPanelMap::Row r(mnid[i],dtcid[i],link[i],station[i],
                         psid[i],plane[i],ppid[i],panel[i],zface[i]);
      ptr->add(r);
    }

    return ptr;
  } // end fromFcl

//-----------------------------------------------------------------------------  
  TrkPanelMapEntity::ptr_t TrkPanelMapMaker::fromDb(TrkPanelMap::cptr_t Table) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl(); 
    int nr = Table->nrow();
    
    for (int i=0; i<nr; i++) {
      ptr->add(Table->rowAt(i));
    }
    
    ptr->print(std::cout);
    
    return ptr;
  }
}
