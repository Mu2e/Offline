#ifndef Mu2eUtilities_resolveTransients_hh
#define Mu2eUtilities_resolveTransients_hh

//
// A utility function to compute values for all of the transients
// in members of a collection.
//
// $Id: resolveTransients.hh,v 1.6 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/Event.h"

namespace mu2e {

  template<typename T>
  void resolveTransients( T const& v, art::Event const& event){

    for ( typename T::const_iterator i=v.begin(), e=v.end();
          i<e; ++i ){
      i->resolveTransients(event);
    }
  }

}

#endif /* Mu2eUtilities_resolveTransients_hh */
