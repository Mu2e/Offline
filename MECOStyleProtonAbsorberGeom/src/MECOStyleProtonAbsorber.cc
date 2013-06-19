//
// Class to represent the system of MECO style Proton Absorber
//
// $Id: MECOStyleProtonAbsorber.cc,v 1.4 2013/06/19 03:41:01 mjlee Exp $
// $Author: mjlee $
// $Date: 2013/06/19 03:41:01 $
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
    _pabs2flag(false),
    //outer PA
    _oPAmaterialName(),
    _oPAzcenter(0),
    _oPAhalflength(0),
    _oPAthickness(0),
    _oPAflag(false)
  {}


  bool MECOStyleProtonAbsorber::isAvailable (int id) const {
    switch (id) {
      case 0:
        return _pabs1flag;
        break;
      case 1:
        return _pabs2flag;
        break;
      case 2:
        return _oPAflag;
      default:
        return false;
    }
    return false;
  }
          

}

