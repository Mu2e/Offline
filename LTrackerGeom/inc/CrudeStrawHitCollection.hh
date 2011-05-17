#ifndef ToyDP_CrudeStrawHitCollection_hh
#define ToyDP_CrudeStrawHitCollection_hh

//
// A collection of CrudeStrawHits, that holds a reference to the
// persistent data and has extra functions:
//   - it provides a view to return a hit by StrawIndex.
//   - it provides a convenience method getStepPointMC().
//
// $Id: CrudeStrawHitCollection.hh,v 1.8 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Rob Kutschke
//
// For all methods for which it makes sense, there is an accessor
// via StrawIndex and an accessor via hit index within the container.
//

// C++ includes.
#include <vector>

// Framework includes.
#include "art/Persistency/Common/Handle.h"

// Mu2e includes.
#include "ToyDP/inc/CrudeStrawHitPData.hh"

namespace art{
  class Event;
}

namespace mu2e {

  class StepPointMC;

  class CrudeStrawHitCollection{

  public:

    // No default constructor by design.
    
    CrudeStrawHitCollection( art::Event const& event,
                             art::Handle<CrudeStrawHitPData> const& hits );

    CrudeStrawHitCollection( art::Event const& event,
                             CrudeStrawHitPData const& hits );

    // Compiler generated versions of the following will be OK:
    //   destructor, copy constructor, assignment operator.

    // For now, this must be called to fix transients lost in persistency.
    void resolveTransients(art::Event const& event );

    // Accessor via index in the persistent container.
    CrudeStrawHit const& get( int i ) const{
      return _hits->at(i);
    }
    CrudeStrawHit const& operator[]( int i ) const{
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

    // Safe version of the previous.
    int strawIndexToHitIndexOrThrow( StrawIndex idx ) const;

  private:

    // Fill the _index variable.
    void FillIndex();

    // These are non-owning pointers.
    art::Event const*         _event;
    CrudeStrawHitPData const* _hits;

    // A second view of the hits, via StrawIndex.
    std::vector<int>  _index;

  };
}

#endif
