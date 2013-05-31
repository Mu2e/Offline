//
// Class to represent the system of MECO style Proton Absorber
//
// $Id: MECOStyleProtonAbsorber.cc,v 1.3 2013/05/31 18:07:18 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:18 $
//
// Original author MyeongJae Lee

#include <iostream>

#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"

namespace mu2e {


  MECOStyleProtonAbsorber::MECOStyleProtonAbsorber() : 
    _parts(),
    _vdHL(0),
    _materialName(),
    _distfromtargetend(0),
    _halflength(0),
    _thickness(0),
    _pabs1flag(false),
    _pabs2flag(false)
  {}


  bool MECOStyleProtonAbsorber::isAvailable (int id) const {
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

