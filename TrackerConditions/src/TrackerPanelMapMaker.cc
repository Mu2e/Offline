// clang-format off
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMapMaker.hh"
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
  TrackerPanelMap::ptr_t TrackerPanelMapMaker::fromFcl() {

    auto ptr = std::make_shared<TrackerPanelMap>();

    std::vector<int> mnid    = config_.mnid   ();
    std::vector<int> dtcid   = config_.dtcid  ();
    std::vector<int> link    = config_.link   ();
    std::vector<int> uniquePlane   = config_.uniquePlane  ();
    std::vector<int> ppid    = config_.ppid   ();
    std::vector<int> panel   = config_.panel  ();
    std::vector<int> zface   = config_.zface  ();

    int npanels = mnid.size();
    for (int i=0; i<npanels; i++) {
      TrkPanelMap::Row r(mnid[i],dtcid[i],link[i],uniquePlane[i],ppid[i],panel[i],zface[i]);
      ptr->add(r);
    }

    return ptr;
  } // end fromFcl

  //-----------------------------------------------------------------------------
  TrackerPanelMap::ptr_t TrackerPanelMapMaker::fromDb(TrkPanelMap::cptr_t Table) {
    auto ptr = std::make_shared<TrackerPanelMap>();
    int nr = Table->nrow();

    for (int i=0; i<nr; i++) {
      ptr->add(Table->rowAt(i));
    }

    ptr->print(std::cout);

    return ptr;
  }
}
