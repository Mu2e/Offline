//
// construct the TrackerStatus object from fcl configuration
//

// C++ includes
#include <iostream>
#include <sstream>
#include <memory>

// Mu2e includes
#include "TrackerConditions/inc/TrackerStatusMaker.hh"
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {

  TrackerStatus::ptr_t TrackerStatusMaker::fromFcl() {

    TrackerStatus::estat_t estats;
    for(auto const& estat : config_.status()) {
      estats.emplace(StrawId(std::get<0>(estat)),
	  StrawIdMask(std::get<1>(estat)), 
	  StrawStatus(std::get<2>(estat)));
    }
    auto const& settings = config_.settings();
    if ( settings.verbose() > 0 ) {
      cout << "TrackerStatus created with " << estats.size() 
	<< " elements " << endl;
      if(settings.verbose() > 1){
	for(auto const& estat : estats) {
	  std::cout << "Tracker Element Id " << estat.sid_
	    << " level " << estat.mask_.levelName()  
	    << " mask " << estat.mask_.mask()  
	    << " status " << estat.status_ << std::endl;
	}
      }
    }
    // make shared_ptr to the TrackerStatus object on the heap
    return std::make_shared<TrackerStatus>(estats);
  }

  TrackerStatus::ptr_t TrackerStatusMaker::fromDb(
      TrkPlaneStatus::cptr_t tpls_p,
      TrkPanelStatus::cptr_t   tpas_p,
      TrkStrawStatusLong::cptr_t   tssl_p,
      TrkStrawStatusShort::cptr_t   tsss_p ) {
    // convert the tables to TrackerElementStatus sets
    TrackerStatus::estat_t estatus;
    auto const& settings = config_.settings();
    if ( settings.verbose() > 1 ) {
      cout << "TrackerStatus fromDb, with tables:" << endl;
      cout << "Plane table has " << tpls_p->rows().size() << " rows " << endl;
      cout << "Panel table has " << tpas_p->rows().size() << " rows " << endl;
      cout << "Long-term Straw table has " << tssl_p->rows().size() << " rows " << endl;
      cout << "Short-term Straw table has " << tsss_p->rows().size() << " rows " << endl;
    }

    for (auto const& row : tpls_p->rows()) estatus.insert(TrackerElementStatus(row.id(),tpls_p->sidMask(),row.status()));
    for (auto const& row : tpas_p->rows()) estatus.insert(TrackerElementStatus(row.id(),tpas_p->sidMask(),row.status()));
    for (auto const& row : tssl_p->rows()) estatus.insert(TrackerElementStatus(row.id(),tssl_p->sidMask(),row.status()));
    for (auto const& row : tsss_p->rows()) estatus.insert(TrackerElementStatus(row.id(),tsss_p->sidMask(),row.status()));

    // check for consistency
    //
    unsigned ntotal = tpls_p->rows().size() + tpas_p->rows().size() + tssl_p->rows().size() + tsss_p->rows().size();
    if(estatus.size() != ntotal){
      throw cet::exception("TrackerStatus BadTable")
	<< "input table size inconsistency "<< ntotal 
	<< " " << estatus.size()  << "\n";

    }

    if ( settings.verbose() > 1 ) {
      cout << "Elements have status " << endl;
      for(auto const& estat : estatus){
	std::cout << "Tracker Element Id " << estat.sid_
	  << " level " << estat.mask_.levelName()  
	  << " mask " << estat.mask_.mask()  
	  << " status " << estat.status_ << std::endl;
      }
    }

    return std::make_shared<TrackerStatus>(estatus);
  }

} // namespace mu2e
