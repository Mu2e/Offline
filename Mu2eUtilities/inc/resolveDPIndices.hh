#ifndef Mu2eUtilities_resolveDPIndices_hh
#define Mu2eUtilities_resolveDPIndices_hh

//
// Several utility functions to resolve a DPIndex, or a
// collection of DPIndex's into pointers to the objects
// that they describe.
//
// $Id: resolveDPIndices.hh,v 1.5 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>

// Framework includes.
#include "art/Framework/Core/Event.h"

// Mu2e includes.
#include "ToyDP/inc/DPIndex.hh"

namespace mu2e {

  // Resolve a single DPIndex.
  template<typename T>
  typename T::value_type const * resolveDPIndex( art::Event const& event,
                                                 DPIndex const&    dpi ){
    
    art::Handle<T> handle;
    event.get( dpi.id, handle);
    
    return &handle->at(dpi.index);
  }

  // Resolve a vector of DPIndices.
  // Assume that ProductIDs may be different from one DPIndex to the other.
  template<typename T>
  void resolveDPIndices( art::Event const&                  event,
                         std::vector<DPIndex> const&        indices,
                         std::vector<typename T::value_type const*>& vout
                         ){

    for ( std::vector<DPIndex>::const_iterator 
            i = indices.begin(),
            e = indices.end();
          i!=e; ++i ){
     
      vout.push_back( resolveDPIndex<T>(event,*i) );
      
    }
    
  }

  // Resolve multiple objects within a single ProductID.
  template<typename T>
  void resolveDPIndices( art::Event const&                  event,
                         art::ProductID const&              id,
                         std::vector<int> const&            indices,
                         std::vector<typename T::value_type const*>& vout
                         ){
    
    art::Handle<T> handle;
    event.get( id, handle);

    for ( std::vector<int>::const_iterator 
            i = indices.begin(),
            e = indices.end();
          i!=e; ++i ){
     
      vout.push_back( &handle->at(*i) );
      
    }
    
  }

}

#endif /* Mu2eUtilities_resolveDPIndices_hh */
