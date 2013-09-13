#ifndef CosmicRayShieldGeom_CRSScintillatorModuleId_hh
#define CosmicRayShieldGeom_CRSScintillatorModuleId_hh

//
// Identifier of one module in CosmicRayShield
//

//
// $Id: CRSScintillatorModuleId.hh,v 1.5 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG; based on Rob Kutschke's SectorId
//

#include <ostream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorShieldId.hh"

namespace mu2e 
{

  class CRSScintillatorModuleId
  {

    friend class CRSScintillatorLayerId;
    friend class CRSScintillatorBarId;

    public:

    CRSScintillatorModuleId():
      _shieldId(CRSScintillatorShieldId(-1)),
      _moduleNumber(-1)
    {
    }
/*
    CRSScintillatorModuleId( CRSScintillatorShieldId shieldId,
                             int moduleNumber):
      _shieldId(shieldId),
      _moduleNumber(moduleNumber)
    {
    }
*/   
    CRSScintillatorModuleId( int shieldNumber,
                             int moduleNumber):
      _shieldId(CRSScintillatorShieldId(shieldNumber)),
      _moduleNumber(moduleNumber)
    {
    }
   

    // Compiler generated d'tor, copy and assignment constructors
    // should be OK.

    const CRSScintillatorShieldId getShieldId() const {return _shieldId;}

    const int getShieldNumber() const {return _shieldId;}
    const int getModuleNumber() const {return _moduleNumber;}

    bool operator==(CRSScintillatorModuleId const & rhs) const
    {
      return ( _shieldId == rhs._shieldId && _moduleNumber == rhs._moduleNumber );
    }

    bool operator!=(CRSScintillatorModuleId const & rhs) const
    {
      return !( *this == rhs);
    }

    private:

    CRSScintillatorShieldId _shieldId;
    int _moduleNumber;
  };

  inline std::ostream& operator<<(std::ostream& ost, const CRSScintillatorModuleId& moduleId )
  {
    ost << moduleId.getShieldId() << " " << moduleId.getModuleNumber();
    return ost;
  }

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorModuleId_hh */
