// clang-format off
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMapMaker.hh"
#include "cetlib_except/exception.h"
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

    size_t npanels = mnid.size();
    if (dtcid.size() != npanels || link.size() != npanels || uniquePlane.size() != npanels ||
        ppid.size()  != npanels || panel.size() != npanels || zface.size() != npanels) {
      throw cet::exception("BADCONFIG")
        << "TrackerPanelMap fcl columns must all have the same length; mnid has " << npanels
        << " entries, dtcid " << dtcid.size() << ", link " << link.size() << ", uniquePlane " << uniquePlane.size()
        << ", ppid " << ppid.size() << ", panel " << panel.size() << ", zface " << zface.size() << "\n";
    }
    for (size_t i=0; i<npanels; i++) {
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
