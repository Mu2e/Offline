#ifndef ToyDP_CrudeStrawHitCollection_hh
#define ToyDP_CrudeStrawHitCollection_hh

//
// A collection of CrudeStrawHits, that holds a reference to the
// persistent data and has extra functions.
//
// $Id: CrudeStrawHitCollection.hh,v 1.1 2009/10/22 15:53:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/22 15:53:23 $
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
    
    CrudeStrawHitCollection( edm::Event const& event,
			     edm::Handle<CrudeStrawHitPData> const& hits );

    CrudeStrawHitCollection( edm::Event const& event,
			     CrudeStrawHitPData const& hits );

    ~CrudeStrawHitCollection();

    // Accessor via index in the presistent container.
    CrudeStrawHit const& get( int i ) const{
      return _hits.at(i);
    }

    // Test if StrawIndex has a hit.
    bool hasHitByStrawIndex( StrawIndex idx ) const{
      return ( _index[idx.asInt()].asInt() != -1 );
    }

    // Accessor via StrawIndex.
    CrudeStrawHit const& getByStrawIndex( StrawIndex idx ) const{
      return _hits.at( _index[idx.asInt()].asInt() );
    }

    // Access the persistent data directly.
    CrudeStrawHitPData const& getPData() const{
      return _hits;
    }

    // Fill the array of pointers to const, v, elements of which point 
    // to the precursors of the specified hit.
    void getStepPointMC( int i,
			 std::vector<StepPointMC const*>& v ) const;

  private:

    // Fill the _index variable.
    void FillIndex();
    
    edm::Event const&         _event;
    CrudeStrawHitPData const& _hits;

    // A second view of the hits, via StrawIndex.
    // These are non-owning pointers.
    std::vector<StrawIndex>   _index;

  };
}

#endif
