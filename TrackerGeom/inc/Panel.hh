#ifndef TrackerGeom_Panel_hh
#define TrackerGeom_Panel_hh
//
// Holds information about one panel in a tracker.
//

//
// $Id: Panel.hh,v 1.14 2013/03/26 23:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/26 23:28:23 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <array>
#include <iomanip>
#include <iostream>

#include "TrackerGeom/inc/Layer.hh"
#include "DataProducts/inc/PanelId.hh"
#include "DataProducts/inc/StrawId2.hh"
#include "GeomPrimitives/inc/TubsParams.hh"


#include "CLHEP/Vector/ThreeVector.h"
#ifndef __CINT__
#include "boost/bind.hpp"
#endif

#include <iostream>

namespace mu2e {

  class Tracker;

  class Panel{


    friend class Plane;
    friend class TTracker;
    friend class TTrackerMaker;

  public:

    Panel():_id(PanelId(-1,-1)),
            _id2(StrawId2())
    {}
    Panel( const PanelId& id ):_id(id){}
    Panel( const PanelId& id, const StrawId2& sid2):_id(id),_id2(sid2){}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const PanelId&  id() const { return _id;}
    const StrawId2& id2() const { return _id2;}

    // const std::vector<Layer>& getLayers() const{
    //   return _layers;
    // }

    const auto& getStrawPointers() const{
      return _straws2_p;
    }

    int nLayers() const{
      return StrawId2::_nlayers;
    }

    int nStraws() const { return StrawId2::_nstraws; }

    const Layer& getLayer ( int n ) const {
      return _layers.at(n);
    }

    const Layer& getLayer ( const LayerId& layid) const {
      return _layers.at(layid.getLayer());
    }

    const Straw& getStraw ( const StrawId& strid ) const{
      return _layers.at(strid.getLayer()).getStraw(strid);
    }

    // this getStraw checks if this is the right panel
    const Straw& getStraw ( const StrawId2& strid2 ) const{
      if ( _id2.samePlane(strid2) && _id2.samePanel(strid2) ) {
        return *(_straws2_p.at((strid2.asUint16() & StrawId2::_strawmsk)));
      } else {
        throw cet::exception("RANGE")
          << __func__ << " Inconsistent straw/panel request " << strid2
          << "/" << id2();
      }
    }

    // Mid-point position of the average (over the layers) of the primary
    // straws, and (collective) straw direction.
    // (The primary straw of each layer is the straw used to establish position.
    //  In the TTracker the primary straw is the innermost straw.)
    // *** In a multi-layer geometry, the straw0MidPoint ***
    // ***        need not lie on any actaul straw       ***
    CLHEP::Hep3Vector straw0MidPoint()  const { return _straw0MidPoint;  }
    CLHEP::Hep3Vector straw0Direction() const { return _straw0Direction; }

    // Formatted string embedding the id of the panel.
    std::string name( std::string const& base ) const;

    // const std::vector<double>& boxHalfLengths() const { return _boxHalfLengths; }

    // double         boxRxAngle()     const { return _boxRxAngle;     }
    // double         boxRyAngle()     const { return _boxRyAngle;     }
    double         boxRzAngle()     const { return _boxRzAngle;     }
    // const CLHEP::Hep3Vector&    boxOffset()      const { return _boxOffset;      }

    std::vector<CLHEP::Hep3Vector> const& getBasePosition() const{
      return _basePosition;
    }

    CLHEP::Hep3Vector const& getBaseDelta() const{
      return _baseDelta;
    }

    TubsParams  const& getEBKeyParams()           const { return _EBKey; }
    TubsParams  const& getEBKeyShieldParams()     const { return _EBKeyShield; }
    std::string const& getEBKeyMaterial()         const { return _EBKeyMaterial; }
    std::string const& getEBKeyShieldMaterial()   const { return _EBKeyShieldMaterial; }
    double             getEBKeyPhiExtraRotation() const { return _EBKeyPhiExtraRotation; }

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker& tracker ) const;

#ifndef __CINT__
    /*
    template <class F>
    void for_each_layer( F f) const{
      std::for_each ( _layers.begin(),
                      _layers.end(),
                      f);
    }

    template <class F>
    void for_each_straw( F f) const {
      for_each_layer( boost::bind( Layer::for_each<F>, _1, f));
    }
    */

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
            i !=e; ++i){
        i->forAllStraws(f);
      }
    }

    template <class F>
    inline void forAllLayers ( F& f) const{
      for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
            i !=e; ++i){
        f(*i);
      }
    }

#endif

  protected:

    PanelId _id;
    StrawId2 _id2;
    std::vector<Layer> _layers;

    std::array<Straw const*, StrawId2::_nstraws> _straws2_p;

    // Vertices of enclosing polygon.
    std::vector<CLHEP::Hep3Vector> corners;

    // Properties of the enclosing logical volume (box).

    // Half lengths of the logical box.
    // std::vector<double> _boxHalfLengths;

    std::vector<CLHEP::Hep3Vector> _basePosition;
    CLHEP::Hep3Vector _baseDelta;

    // Rotations and offsets to place the logical box.
    // placedshape = ( offset + RZ*RX*RY*shape );
    //
    // double _boxRxAngle;
    // double _boxRyAngle;
    double _boxRzAngle; // is it really used? needed?
    // CLHEP::Hep3Vector _boxOffset;

    // Position (in tracker coordinates) of the midpoint, and direction
    // of the average of the primary straw.  Mutable because these are set
    // by fillPointers.
    // TODO -- there is clearly a way to design this such that these mutable
    //         declarations can go away.
    mutable CLHEP::Hep3Vector _straw0MidPoint;
    mutable CLHEP::Hep3Vector _straw0Direction;

    // electronic boards; they are placed wrt to the panel, but in a
    // different mother volume one per panel

    TubsParams  _EBKey;
    std::string _EBKeyMaterial;
    TubsParams  _EBKeyShield;
    std::string _EBKeyShieldMaterial;
    double      _EBKeyPhiExtraRotation;

  };

}  //namespace mu2e
#endif /* TrackerGeom_Panel_hh */
