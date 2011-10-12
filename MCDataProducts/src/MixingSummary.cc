//
// Status information from one call to an event mixing module.
//
// $Id: MixingSummary.cc,v 1.1 2011/10/12 20:01:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/10/12 20:01:17 $
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/MixingSummary.hh"

#include <iostream>

using namespace std;

namespace mu2e{

  void MixingSummary::print ( std::ostream& ost ) const{
    ost << "Mixing Summary: Overall status: "
        << status_ << endl;

    ost << "Number of event Ids: " << eventIDs_.size() << "\n";

    for ( art::EventIDSequence::const_iterator i0=eventIDs_.begin(), i=i0, e=eventIDs_.end();
          i != e; ++i ){
      size_t j=i-i0;
      ost << "   EventId: "
          << *i << " Sizes/Deltas: "
          << genSizes_.at(j) << " "
          << simDeltas_.at(j) << " ";
      for ( size_t k=0; k<stepSizes_.size(); ++k){
        if ( j < stepSizes_.at(k).size() ) {
          ost << stepSizes_.at(k).at(j) << " ";
        } else{
          if ( k != StepInstanceName::unknown ) {
            ost << "-- ";
          }
        }

      }
      ost << " | " << pointTrajectoryDeltas_.at(j) << "\n";
      ost << "            Status: "
          << eventStatus_.at(j)
          << endl;
    }

  } // end MixingSummary::print

} // end namespace mu2e
