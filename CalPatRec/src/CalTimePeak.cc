//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  Original author Pavel Murat
//

#include "CalPatRec/inc/CalTimePeak.hh"

using namespace std;

namespace mu2e {

  // Constructors
  CalTimePeak::CalTimePeak() {
    _cluster  = NULL;
    _cprIndex = -1;
    _x        = 0.;
    _y        = 0.;
    _z        = 0.;
    _tpeak    = 0.;
    _shcol    = NULL;
    _shfcol   = NULL;
  }

  CalTimePeak::CalTimePeak(const CaloCluster* Cl, double X, double Y, double Z) {
    _cluster  = Cl;
    _x        = X;
    _y        = Y;
    _z        = Z;
    _cprIndex = -1;
    _tpeak    = 0;
    _shcol    = NULL;
    _shfcol   = NULL;
  }

  CalTimePeak::~CalTimePeak() {
  }

  void CalTimePeak::clear () {
    _cluster  = NULL;
    _cprIndex = -1;
    _x        = 0.;
    _y        = 0.;
    _z        = 0;
    _tpeak    = 0;
    _shcol    = NULL;
    _shfcol   = NULL;
    _index.clear();
  }

} // end namespace mu2e


