#ifndef Mu2eUtilities_resolveDPIndices_HH
#define Mu2eUtilities_resolveDPIndices_HH

//
// Several utility functions to resolve a DPIndex, or a
// collection of DPIndex's into pointers to the objects
// that they describe.
//
// $Id: resolveDPIndices.hh,v 1.1 2009/10/22 15:52:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/22 15:52:07 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"

// Mu2e includes.
#include "ToyDP/inc/DPIndex.hh"

namespace mu2e {

  template<typename T>
  typename T::value_type const * resolveDPIndex( edm::Event const& event,
						 DPIndex const&    dpi ){
    
    edm::Handle<T> handle;
    event.get( dpi.id, handle);
    
    return &handle->at(dpi.index);
  }

  template<typename T>
  void resolveDPIndices( edm::Event const&                  event,
			 std::vector<DPIndex> const&        indices,
			 std::vector<typename T::value_type const*>& vout
			 ){

    for ( std::vector<DPIndex>::const_iterator 
	    i = indices.begin(),
	    e = indices.end();
	  i!=e; ++i ){
     
      // Get a handle to a collection described by a ProductID.
      // edm::Handle<T> handle;
      //event.get( i->id, handle);
      
      //vout.push_back( &handle->at(i->index) );
      //DPIndex const& idx(*i);
      vout.push_back( resolveDPIndex<T>(event,*i) );
      
    }
    
  }

}

#endif
