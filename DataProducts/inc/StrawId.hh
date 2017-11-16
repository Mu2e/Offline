#ifndef TrackerGeom_StrawId_hh
#define TrackerGeom_StrawId_hh
//
// Identifier of one straw in a tracker.
//

//
// $Id: StrawId.hh,v 1.7 2011/05/20 19:18:45 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 19:18:45 $
//
// Original author Rob Kutschke
//
#include <ostream>
#include <iostream>
#include "DataProducts/inc/LayerId.hh"
#include <string>
// art includes
#include "cetlib_except/exception.h"

namespace mu2e {

  class StrawId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const StrawId& s );

class StrawId{

public:

  StrawId(std::string const& asstring);

  StrawId():
    _lid(),
    _n(-1){
  }

  StrawId( LayerId layer,
           int n
           );

  StrawId( PanelId panelid,
           int layer,
           int n
           ):
    _lid(panelid,layer),
    _n(n){

    if ( n%2!=layer%2 ) {
      std::cerr << "CONFIG " 
        //      throw cet::exception("CONFIG")
                << "StrawId(PanelId, int, int): incorrect straw in layer "
                << panelid << "_" << layer << "_" << n
                << "\n";
    }

  }

  StrawId( PlaneId plane,
           int panel,
           int layer,
           int n
           ):
    _lid(LayerId(plane,panel,layer)),
    _n(n){

    if ( n%2!=layer%2 ) {
      std::cerr << "CONFIG "
        //      throw cet::exception("CONFIG")
                << "StrawId(PlaneId, int, int, int): incorrect straw in layer "
                << plane << "_" << panel << "_" << layer <<  "_" << n
                << "\n";
    }
    
  }

  // Use compiler-generated copy c'tor, copy assignment, and d'tor.

  const PlaneId& getPlaneId() const {
    return _lid.getPlaneId();
  }

  const PanelId& getPanelId() const {
    return _lid.getPanelId();
  }

  const LayerId& getLayerId() const {
    return _lid;
  }

  int getPlane() const{
    return _lid.getPlane();
  }

  int getPanel() const{
    return _lid.getPanel();
  }

  int getLayer() const{
    return _lid.getLayer();
  }

  int getStraw() const{
    return _n;
  }

  int getStation() const{
    return _lid.getStation();
  }

  bool operator==( StrawId const& rhs) const{
    return ( _lid == rhs._lid && _n == rhs._n );
  }

  bool operator!=( StrawId const& rhs) const{
    return !( *this == rhs);
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const StrawId& s ){
    ost << s.getLayerId() << "_"
        << s._n;
    return ost;
  }

private:

  LayerId _lid;
  int     _n;

};

}
#endif /* TrackerGeom_StrawId_hh */
