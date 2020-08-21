//
// Representation of one ScintillatorShield in CosmicRayShield
//
//
// Original author KLG; somewhat based on Rob Kutschke's Device
//

#include <string>
#include <vector>
#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

using namespace std;

namespace mu2e {

  CRSScintillatorShield::CRSScintillatorShield(CRSScintillatorShieldId const & id,
                                               std::string const & name,
                                               const std::shared_ptr<CRSScintillatorBarDetail> barDetails,
                                               const std::string &absorberMaterialName, const std::string &aluminumSheetMaterialName, const std::string &FEBMaterialName,
                                               CRSScintillatorShieldId precedingSector, int sectorType, int countersPerModule) :
    _id(id),
    _name(name),
    _barDetails(barDetails),
    _absorberMaterialName(absorberMaterialName),
    _aluminumSheetMaterialName(aluminumSheetMaterialName),
    _FEBMaterialName(FEBMaterialName),
    _precedingSector(precedingSector),
    _sectorType(sectorType),
    _countersPerModule(countersPerModule)
  {
  }

  string CRSScintillatorShield::name( string const& base ) const
  {
    ostringstream os;

    os << base
       << _id;
    return os.str();
  }

}
