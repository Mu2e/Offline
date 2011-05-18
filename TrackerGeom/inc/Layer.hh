#ifndef TrackerGeom_Layer_hh
#define TrackerGeom_Layer_hh
//
// Hold information about one Layer in a tracker.
//
// $Id: Layer.hh,v 1.6 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <deque>

#include "TrackerGeom/inc/LayerId.hh"
#include "TrackerGeom/inc/Straw.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Tracker;

  class Layer{

    friend class Sector;
    friend class Device;
    friend class LTracker;
    friend class TTracker;
    friend class LTrackerMaker;
    friend class TTrackerMaker;

  public:

    // A free function, returning void, that takes a const Layer& as an argument.
    typedef void (*LayerFunction)( const Layer& s);

    Layer();

    Layer(const LayerId& id);

    Layer(const LayerId&   id,
          int        nStraws,
          const CLHEP::Hep3Vector& origin,
          const CLHEP::Hep3Vector& delta
          );

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const LayerId& Id() const { return _id;}

    int nStraws() const { return _straws.size(); }

    const CLHEP::Hep3Vector getOrigin() const {return _orig;}
    const CLHEP::Hep3Vector getDelta()  const {return _delta;}

    const Straw& getStraw( int n ) const {
      return *_straws.at(n);
    }

    const Straw& getStraw( const StrawId& id ) const {
      return getStraw(id._n);
    }

    const std::vector<const Straw*>& getStraws() const {
      return _straws;
    }

    // Formatted string embedding the id of the layer.
    std::string name( std::string const& base ) const;

    // Return Id of the last straw in the layer.
    // Return an illegal id if there are no straws.
    //  - which should never happen.
    StrawId getIdLastStraw() const{

      return ( _straws.size() != 0 )?
        _straws.back()->Id():
        StrawId();
    }

    // Options:
    // 1) return F or void; std-like says return F.
    //    Can be inefficient if F has a lot of state.
    //    THis is a copy in copy out model.
    // 2) Order is not defined.  In practice when run on
    //    vector it is defined.
    template <class F>
    void for_each( F f ) const {
      std::for_each( _straws.begin(),
                     _straws.end(),
                     f);
    }

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( std::vector<const Straw*>::const_iterator i=_straws.begin(), e=_straws.end();
            i !=e; ++i){
        f(**i);
      }
    }

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

  protected:

    LayerId _id;

    // Number of straws.  Needed because of 2 phase construction.
    // The member _straws is not filled until the second phase
    // but this is neede beforehand. Keep it strictly private.
    int32_t _nStraws;

    // Nominal position of wire 0 and offset from wire 0 to wire 1.
    // This is exactly only all wires are in their nominal positions.
    CLHEP::Hep3Vector _orig;
    CLHEP::Hep3Vector _delta;

    // Pointers to the straws in this layer.
    // These pointers do not own the straws to which they point.
    // These are not persisted and may need to be recomputed after readback.
    mutable std::vector<const Straw*> _straws;

    std::vector<StrawIndex> _indices;

  };
}

#endif /* TrackerGeom_Layer_hh */
