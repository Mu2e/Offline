#ifndef TrackerGeom_Plane_hh
#define TrackerGeom_Plane_hh
//
// Hold information about one plane in a tracker.
//
// $Id: Plane.hh,v 1.10 2014/04/11 04:39:13 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:39:13 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>

// Mu2e includes
#include "DataProducts/inc/PlaneId.hh"
#include "TrackerGeom/inc/Panel.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Tracker;

  class Plane{

    friend class TTracker;
    friend class TTrackerMaker;

  public:

    // A free function, returning void, that takes a const Plane& as an argument.
    typedef void (*PlaneFunction)( const Plane& s);

    Plane():_id(PlaneId()),_rotation(0.),_origin(),_panels(),_exists(true){}

    explicit Plane( const PlaneId& id,
            CLHEP::Hep3Vector const& origin = CLHEP::Hep3Vector(0.,0.,0.),
            double rotation = 0., bool exists = true ):
      _id(id),
      _rotation(rotation),
      _origin(origin),
      _exists(exists) {
    }

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const PlaneId&  id()  const { return _id;}

    double rotation() const { return _rotation; }

    const CLHEP::Hep3Vector & origin() const { return _origin; }

    int nPanels() const{
      return _panels.size();
    }

    const std::array<Panel,StrawId::_npanels>& getPanels () const{
      return _panels;
    }

    const Panel& getPanel ( int n) const {
      return _panels.at(n);
    }

    const Panel& getPanel ( const PanelId& pnlid ) const{
      return _panels.at(pnlid.getPanel());
    }

    const Straw& getStraw ( const StrawId& strid ) const{
      return _panels.at(strid.getPanel()).getStraw(strid);
    }

    // Formatted string embedding the id of the panel.
    std::string name( std::string const& base ) const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

    bool exists() const {
      return _exists;
    }

#ifndef __CINT__

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( const auto& panel : _panels ){
        panel.forAllStraws(f);
      }
    }

    template <class F>
    inline void forAllPanels ( F& f) const{
      for ( const auto& panel : _panels ) {
        f(panel);
      }
    }

#endif

  protected:

    PlaneId             _id;
    double              _rotation;
    CLHEP::Hep3Vector   _origin;
    std::array<Panel,StrawId::_npanels> _panels;
    bool                _exists;
  };

} //namespace mu2e

#endif /* TrackerGeom_Plane_hh */
