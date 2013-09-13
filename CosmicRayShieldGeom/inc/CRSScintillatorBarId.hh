#ifndef CosmicRayShieldGeom_CRSScintillatorBarId_hh
#define CosmicRayShieldGeom_CRSScintillatorBarId_hh
//
// Identifier of one CRSScintillatorBar in CosmicRayShield.
//

//
// $Id: CRSScintillatorBarId.hh,v 1.5 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG somewhat based on Rob Kutschke's StrawId
//
#include <ostream>
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayerId.hh"

namespace mu2e 
{

  class CRSScintillatorBarId
  {

    public:

    CRSScintillatorBarId():
      _layerId(CRSScintillatorLayerId()),
      _barNumber(-1)
    {
    }

    CRSScintillatorBarId( CRSScintillatorLayerId layerId,
                          int barNumber):
      _layerId(layerId),
      _barNumber(barNumber)
    {
    }

    CRSScintillatorBarId( int shieldNumber,
                          int moduleNumber,
                          int layerNumber,
                          int barNumber):
      _layerId(CRSScintillatorLayerId(shieldNumber,moduleNumber,layerNumber)),
      _barNumber(barNumber)
    {
    }

    // Compiler generated d'tor, copy c'tor and assignment
    // operators should be should be OK.

    const CRSScintillatorShieldId& getShieldId() const {return _layerId._moduleId._shieldId;}
    const CRSScintillatorModuleId& getModuleId() const {return _layerId._moduleId;}
    const CRSScintillatorLayerId& getLayerId() const {return _layerId;}

    const int getShieldNumber() const {return _layerId._moduleId._shieldId;}
    const int getModuleNumber() const {return _layerId._moduleId._moduleNumber;}
    const int getLayerNumber() const {return _layerId._layerNumber;}
    const int getBarNumber() const {return _barNumber;}

    bool operator==( CRSScintillatorBarId const & rhs) const 
    {
      return ( _layerId == rhs._layerId && _barNumber == rhs._barNumber );
    }

    bool operator!=( CRSScintillatorBarId const & rhs) const 
    {
      return !( *this == rhs);
    }

    private:

    CRSScintillatorLayerId _layerId;
    int _barNumber;

  };

  inline std::ostream& operator<<(std::ostream& ost, const CRSScintillatorBarId& barId)
  {
    ost << "CRSScintillatorBar Id: ("
        << barId.getLayerId() << " "
        << barId.getBarNumber()
        << " )";
    return ost;
  }

}
#endif /* CosmicRayShieldGeom_CRSScintillatorBarId_hh */
