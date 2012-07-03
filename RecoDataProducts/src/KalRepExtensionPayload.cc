//
// The extended persistent information from a KalRep.
//
// $Id: KalRepExtensionPayload.cc,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Contact person Rob Kutschke
//

#include "RecoDataProducts/inc/KalRepExtensionPayload.hh"

#include <iostream>

namespace mu2e {

  void KalRepExtensionPayload::print( std::ostream& ost, bool doEndl ) const{
    if ( doEndl ) ost << std::endl;
  }

}  // namespace mu2e
