#ifndef CRSScintillatorBarId_HH
#define CRSScintillatorBarId_HH
//
// Identifier of one CRSScintillatorBar in CosmicRayShield.
//

//
// $Id: CRSScintillatorBarId.hh,v 1.1 2011/03/09 19:22:03 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:22:03 $
//
// Original author KLG somewhat based on Rob Kutschke's StrawId
//
#include <ostream>
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayerId.hh"

namespace mu2e { 

  class CRSScintillatorBarId{

  public:

    CRSScintillatorBarId():
      _layerId(CRSScintillatorLayerId()),
      _barNumber(-1){
    }
  
    CRSScintillatorBarId( CRSScintillatorLayerId layerId,
                          int barNumber
                          ):
      _layerId(layerId),
      _barNumber(barNumber){
    }
  
    CRSScintillatorBarId( CRSScintillatorModuleId moduleId,
                          int layerNumber,
                          int barNumber
                          ):
      _layerId(moduleId,layerNumber),
      _barNumber(barNumber){
    }

    CRSScintillatorBarId( CRSScintillatorShieldId shieldId,
                          int moduleNumber,
                          int layerNumber,
                          int barNumber
                          ):
      _layerId(CRSScintillatorLayerId(shieldId,moduleNumber,layerNumber)),
      _barNumber(barNumber){
    }

    // Compiler generated d'tor, copy c'tor and assignment
    // operators should be should be OK.

    const CRSScintillatorShieldId& getShieldId() const {
      return _layerId._moduleId._shieldId;
    }

    const CRSScintillatorModuleId& getModuleId() const {
      return _layerId._moduleId;
    }

    const CRSScintillatorLayerId& getLayerId() const {
      return _layerId;
    }
  
    const int getShieldNumber() const {
      return _layerId._moduleId._shieldId; // it is a typdef for now anyway
    }

    const int getModuleNumber() const {
      return _layerId._moduleId._moduleNumber;
    }

    const int getLayerNumber() const {
      return _layerId._layerNumber;
    }

    const int getBarNumber() const {
      return _barNumber;
    }

    bool operator==( CRSScintillatorBarId const & rhs) const {
      return ( _layerId == rhs._layerId && _barNumber == rhs._barNumber );
    }

    bool operator!=( CRSScintillatorBarId const & rhs) const {
      return !( *this == rhs);
    }

  private:

    CRSScintillatorLayerId _layerId;
    int32_t _barNumber;
  
  };

  inline std::ostream& operator<<(std::ostream& ost, 
                                  const CRSScintillatorBarId& barId ){
    ost << "CRSScintillatorBar Id: ("
        << barId.getLayerId() << " "
        << barId.getBarNumber()
        << " )";
    return ost;
  }

}
#endif
