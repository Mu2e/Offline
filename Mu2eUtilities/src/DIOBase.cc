// DIOBase base function 
//
// $Id: DIOBase.cc,v 1.1 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

// Framework includes
#include "cetlib/exception.h"

namespace mu2e {

  void DIOBase::checkTable() const {
    
    double valueToCompare = (_table.at(0).first) + 1e9; 
    //order check
    for ( const auto & it : _table ) {
      if (it.first >= valueToCompare) {
        throw cet::exception("Format")
          << "Wrong value in table: " << it.first;
      }
    }
  }

}

