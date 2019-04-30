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

#ifndef __CINT__
#include "boost/bind.hpp"
#endif

#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

#include "TrackerGeom/inc/Straw.hh"
#include "DataProducts/inc/PanelId.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

namespace mu2e {

  class Panel{

    friend class TrackerMaker;
    friend class Tracker; // needed for deep copy

  public:

    Panel():_id(PanelId()),_straws2_p(){}
    Panel( const PanelId& id ):_id(id),_straws2_p(){}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const PanelId& id() const { return _id;}

    // const std::vector<Layer>& getLayers() const{
    //   return _layers;
    // }

    const auto& getStrawPointers() const{
      return _straws2_p;
    }

    int nLayers() const{
      return StrawId::_nlayers;
    }

    int nStraws() const { return StrawId::_nstraws; }

    const Straw& getStraw( uint16_t n ) const {
      return *(_straws2_p.at(n));
    }

    // const Straw& getStraw ( const StrawId& strid ) const{
    //   return _layers.at(strid.getLayer()).getStraw(strid);
    // }

    // this getStraw checks if this is the right panel
    const Straw& getStraw ( const StrawId& strid2 ) const{
      if ( _id.samePlane(strid2) && _id.samePanel(strid2) ) {
        return *(_straws2_p.at((strid2.asUint16() & StrawId::_strawmsk)));
      } else {
        std::ostringstream msg;
        msg << __func__ << " Inconsistent straw/panel request "
            << strid2 << " / " << _id << std::endl;
        throw cet::exception("RANGE") << msg.str();
      }
    }

    // Mid-point position of the average (over the layers) of the primary
    // straws, and (collective) straw direction.
    // (The primary straw of each layer is the straw used to establish position.
    //  In the Tracker the primary straw is the innermost straw.)
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
    void fillPointers ( const Tracker* tracker ) const;


  protected:

    PanelId _id;

    std::array<Straw const*, StrawId::_nstraws> _straws2_p;

    // Vertices of enclosing polygon.
    std::vector<CLHEP::Hep3Vector> corners;

    // Properties of the enclosing logical volume (box).

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
