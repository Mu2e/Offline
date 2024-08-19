//
// Convert coordinates between DUSAF and Mu2e coordinate systems
// See Mu2e-doc-4741-v6
//

#include "Offline/GeometryService/inc/DUSAFMu2eConverter.hh"

namespace mu2e {

  DUSAFMu2eConverter::DUSAFMu2eConverter():
    _tieIn_inDUSAF_usft( 99242.290578, 99026.647809, 729.758494 ),   // From Table 4
    _tieIn_inMu2e_mm(5794.738, 345.131, 1600.553 ),                  // From Table 4
    _toDUSAF(CLHEP::HepRotation().setRows( {0.54671783, -0.00003629,  0.83731691},    // From equation 11
                                           {0.83731691, -0.00004662, -0.54671784},
                                           {0.00005887,  1.00000000,  0.00000490} )),
    _toMu2e(_toDUSAF.inverse())
  {
  }

  // Equation 11
  CLHEP::Hep3Vector DUSAFMu2eConverter::Mu2e_to_DUSAF( CLHEP::Hep3Vector const & v ) const{
    auto result = _tieIn_inDUSAF_usft  + _toDUSAF*(v-_tieIn_inMu2e_mm)/us_foot_mm;;
    return result;
  }

  // Equation 10
  CLHEP::Hep3Vector DUSAFMu2eConverter::DUSAF_to_Mu2e( CLHEP::Hep3Vector const & v ) const{
    auto result = _tieIn_inMu2e_mm + _toMu2e*(v-_tieIn_inDUSAF_usft)*us_foot_mm;
    return result;
  }

}
