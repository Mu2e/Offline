//
// construct the TrackerStatus object from fcl configuration
//

// C++ includes
#include <iostream>
#include <sstream>
#include <memory>

// Mu2e includes
#include "Offline/TrackerConditions/inc/TrackerStatusMaker.hh"
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {

  TrackerStatus::ptr_t TrackerStatusMaker::fromFcl() {
    auto trkstatptr = std::make_shared<TrackerStatus>();

    for(auto const& estat : config_.status()) {
      trkstatptr->addStatus(StrawId(std::get<0>(estat)),
          StrawIdMask(std::get<1>(estat)),
          StrawStatus(std::get<2>(estat)));
    }
    auto const& settings = config_.settings();
    if ( settings.verbose() > 0 ) {
      cout << "TrackerStatus create from fcl with " << config_.status().size()
        << " elements " << endl;
    }
    if(settings.verbose() > 1){
      trkstatptr->print(cout);
    }
    return trkstatptr;
  }

  TrackerStatus::ptr_t TrackerStatusMaker::fromDb(
      TrkPlaneStatus::cptr_t tpls_p,
      TrkPanelStatus::cptr_t   tpas_p,
      TrkStrawStatusLong::cptr_t   tssl_p,
      TrkStrawStatusShort::cptr_t   tsss_p ) {
    // create return object
    auto trkstatptr = std::make_shared<TrackerStatus>();
    auto const& settings = config_.settings();
    if ( settings.verbose() > 1 ) {
      cout << "TrackerStatus fromDb, with tables:" << endl;
      cout << "Plane table has " << tpls_p->rows().size() << " rows " << endl;
      cout << "Panel table has " << tpas_p->rows().size() << " rows " << endl;
      cout << "Long-term Straw table has " << tssl_p->rows().size() << " rows " << endl;
      cout << "Short-term Straw table has " << tsss_p->rows().size() << " rows " << endl;
    }

    for (auto const& row : tpls_p->rows()) trkstatptr->addStatus(row.id(),tpls_p->sidMask(),row.status());
    for (auto const& row : tpas_p->rows()) trkstatptr->addStatus(row.id(),tpas_p->sidMask(),row.status());
    for (auto const& row : tssl_p->rows()) trkstatptr->addStatus(row.id(),tssl_p->sidMask(),row.status());
    for (auto const& row : tsss_p->rows()) trkstatptr->addStatus(row.id(),tsss_p->sidMask(),row.status());

    if ( settings.verbose() > 1 ) trkstatptr->print(cout);
    return trkstatptr;
  }

} // namespace mu2e
