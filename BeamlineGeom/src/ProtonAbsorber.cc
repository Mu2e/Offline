//
// Class to represent the system of MECO style Proton Absorber
//
//
// Original author MyeongJae Lee

#include <iostream>

#include "Offline/BeamlineGeom/inc/ProtonAbsorber.hh"

namespace mu2e {


  ProtonAbsorber::ProtonAbsorber() :
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
    _oPA1flag(false),
    _oPA2flag(false)
  {}


  bool ProtonAbsorber::isAvailable (int id) const {
    switch (id) {
      case ProtonAbsorberId::pabs1 :
        return _pabs1flag;
        break;
      case ProtonAbsorberId::pabs2 :
        return _pabs2flag;
        break;
      case ProtonAbsorberId::opabs1 :
        return _oPA1flag;
      case ProtonAbsorberId::opabs2 :
        return _oPA2flag;
      default:
        return false;
    }
    return false;
  }


}
