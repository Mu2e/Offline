#ifndef CosmicRayShieldGeom_CRSScintillatorLayerId_hh
#define CosmicRayShieldGeom_CRSScintillatorLayerId_hh

//
// Identifier of one Scintillator Layer in CosmicRayShield.
//

//
// $Id: CRSScintillatorLayerId.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author KLG; based on Rob Kutschke LayerId
//

#include <ostream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorModuleId.hh"

namespace mu2e {

  class CRSScintillatorLayerId{

    friend class CRSScintillatorBarId;

  public:

    CRSScintillatorLayerId():
      _moduleId(CRSScintillatorModuleId()),
      _layerNumber(-1){
    }

    CRSScintillatorLayerId( CRSScintillatorModuleId moduleId,
                            int layerNumber
                            ):
      _moduleId(moduleId),
      _layerNumber(layerNumber){
    }

    CRSScintillatorLayerId( CRSScintillatorShieldId shieldId,
                            int moduleNumber,
                            int layerNumber
                            ):
      _moduleId(CRSScintillatorModuleId(shieldId,moduleNumber)),
      _layerNumber(layerNumber){
    }

    // Compiler generated d'tor, copy and assignment constructors
    // should be OK.

    const CRSScintillatorShieldId getShieldId() const{
      return _moduleId._shieldId;
    }
    const CRSScintillatorModuleId getModuleId() const{
      return _moduleId;
    }

    const int getShieldNumber() const{
      return _moduleId._shieldId; // it is a typdef for now anyway
    }

    const int getModuleNumber() const{
      return _moduleId._moduleNumber;
    }

    const int getLayerNumber() const{
      return _layerNumber;
    }

    bool operator==(CRSScintillatorLayerId const & rhs) const{
      return ( _moduleId == rhs._moduleId && _layerNumber == rhs._layerNumber );
    }

    bool operator!=(CRSScintillatorLayerId const & rhs) const{
      return !( *this == rhs);
    }
  
  private:

    CRSScintillatorModuleId _moduleId;
    int32_t _layerNumber;
  
  };

  inline std::ostream& operator<<(std::ostream& ost, 
                                  const CRSScintillatorLayerId& layerId ){
    ost << layerId.getModuleId() << " " << layerId.getLayerNumber();
    return ost;
  }

} //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorLayerId_hh */
