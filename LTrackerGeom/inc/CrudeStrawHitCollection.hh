#ifndef ToyDP_CrudeStrawHitCollection_hh
#define ToyDP_CrudeStrawHitCollection_hh

//
// A collection of CrudeStrawHits, that holds a reference to the
// persistent data and has extra functions:
//   - it provides a view to return a hit by StrawIndex.
//   - it provides a convenience method getStepPointMC().
//
// $Id: CrudeStrawHitCollection.hh,v 1.3 2009/10/28 13:36:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/28 13:36:50 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <vector>

// Framework includes.
#include "DataFormats/Common/interface/Handle.h"

// Mu2e includes.
#include "LTrackerGeom/inc/StrawIndex.hh"
#include "ToyDP/inc/CrudeStrawHitPData.hh"

namespace edm{
  class Event;
}

namespace mu2e {

  class StepPointMC;

  class CrudeStrawHitCollection{

  public:

    // No default constructor by design.
    
    CrudeStrawHitCollection( edm::Event const& event,
			     edm::Handle<CrudeStrawHitPData> const& hits );

    CrudeStrawHitCollection( edm::Event const& event,
			     CrudeStrawHitPData const& hits );

    // Compiler generated versions of the following will be OK:
    //   destructor, copy constructor, assignment operator.

    // Accessor via index in the persistent container.
    CrudeStrawHit const& get( int i ) const{
      return _hits->at(i);
    }

    // Test if StrawIndex has a hit.
    bool hasHitByStrawIndex( StrawIndex idx ) const{
      return ( _index[idx.asInt()] != -1 );
    }

    // Accessor via StrawIndex.
    CrudeStrawHit const& getByStrawIndex( StrawIndex idx ) const{
      return _hits->at( _index[idx.asInt()] );
    }

    // Return hit index in this container, addressed by straw index.
    int indexByStrawIndex( StrawIndex idx ) const{
      return _index[idx.asInt()];
    }

    // Access the persistent data directly.
    CrudeStrawHitPData const& getPData() const{
      return *_hits;
    }

    // Fill the array of pointers to const, v, elements of which point 
    // to the precursors of the specified hit.
    void getStepPointMC( int i,
			 std::vector<StepPointMC const*>& v ) const;

  private:

    // Fill the _index variable.
    void FillIndex();

    // These are non-owning pointers.
    edm::Event const*         _event;
    CrudeStrawHitPData const* _hits;

    // A second view of the hits, via StrawIndex.
    std::vector<int>  _index;

  };
}

#endif
