//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
// Framework includes.
#include "cetlib_except/exception.h"
// C++ includes
#include <ostream>

using namespace std;

namespace mu2e {

  // Default constructor is required for persistable classes
  StrawDigiMC::StrawDigiMC()
    : _strawid(StrawId::_invalid)
    , _validity(Invalid)
  {}

  StrawDigiMC::StrawDigiMC(StrawId sid, PA cpos, FA ctime, FA wetime, SGSPA sgs, Validity validity):
    _strawid(sid), _cpos(cpos), _ctime(ctime), _wtime(wetime), _sgspa(sgs)
    , _validity(validity)
  {}

  StrawDigiMC::StrawDigiMC(const StrawDigiMC& rhs, SGSPA sgspa ) : StrawDigiMC(rhs)  {
    _sgspa = sgspa;
  }

  StrawDigiMC::StrawDigiMC(const StrawDigiMC& rhs, Validity validity): StrawDigiMC(rhs){
    _validity = validity;
  }

  bool StrawDigiMC::isCrossTalk(StrawEnd strawend) const {
    bool retval(false);
    if(!_sgspa[strawend].isNull()){
      retval = _strawid != _sgspa[strawend]->strawId();
    }
    return retval;
  }

  double StrawDigiMC::energySum() const {
    if(_sgspa[0] == _sgspa[1] )
      return _sgspa[0]->ionizingEdep();
    else
      return _sgspa[0]->ionizingEdep() + _sgspa[1]->ionizingEdep();
  }

  double StrawDigiMC::triggerEnergySum(StrawEnd strawend) const {
    return _sgspa[strawend]->ionizingEdep();
  }

  // Print the information found in this hit.
  void StrawDigiMC::print( ostream& ost, bool doEndl ) const {
    ost << "Straw Digi MC Truth for straw ends " << StrawEnd(StrawEnd::cal) << " : " << StrawEnd(StrawEnd::hv)
      << " cluster times : "      << _ctime[0] << " : " << _ctime[1]
      << " Energy: " << energySum();

    if ( doEndl ){
      ost << endl;
    }
  }
}
