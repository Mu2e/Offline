//
// Class to represent the system of MECO style Proton Absorber
//
// $Id: MECOStyleProtonAbsorber.cc,v 1.1 2012/05/10 23:40:59 mjlee Exp $
// $Author: mjlee $
// $Date: 2012/05/10 23:40:59 $
//
// Original author MyeongJae Lee

#include <iostream>

#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"

namespace mu2e {


  MECOStyleProtonAbsorber::MECOStyleProtonAbsorber() : 
    _parts(),
    _ds2zcenter(0),
    _vdHL(0),
    _materialName(),
    _distfromtargetend(0),
    _halflength(0),
    _thickness(0),
    _pabs1flag(false),
    _pabs2flag(false)
  {}


  bool MECOStyleProtonAbsorber::isAvailable (int id) {
    switch (id) {
      case 0:
        return _pabs1flag;
        break;
      case 1:
        return _pabs2flag;
        break;
      default:
        return false;
    }
    return false;
  }
          

}

