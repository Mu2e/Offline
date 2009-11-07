#ifndef Mu2eUtilities_resolveTransients_HH
#define Mu2eUtilities_resolveTransients_HH

//
// A utility function to clear or restore all of the transients
// in a collection.
//
// $Id: resolveTransients.hh,v 1.1 2009/11/07 01:12:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/07 01:12:57 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "FWCore/Framework/interface/Event.h"

namespace mu2e {

  // For all elements in a container, resolve their indices.
  template<typename T>
  void resolveTransients( T const& v, edm::Event const& event){
    
    for ( typename T::const_iterator i=v.begin(), e=v.end();
	  i<e; ++i ){
      
      // This is needed until we move to a more modern version
      // of genreflex that knows how to set transients on readback.
      i->resetTransients();

      // This is the real work.
      i->resolveTransients(event);
    }
  }

  // For all elements in a container, reset their transients to
  // the not-yet-defined state.  Eventually the edm will do this auto-magically.
  template<typename T>
  void resetTransients( T const& v ){
    
    for ( typename T::const_iterator i=v.begin(), e=v.end();
	  i<e; ++i ){
      
      // This is needed until we move to a more modern version
      // of genreflex that knows how to set transients on readback.
      i->resetTransients();

    }
  }


}

#endif
