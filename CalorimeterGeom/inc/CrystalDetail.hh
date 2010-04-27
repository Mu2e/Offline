#ifndef CRYSTALDETAIL_HH
#define CRYSTALDETAIL_HH
// $Id: CrystalDetail.hh,v 1.3 2010/04/27 18:46:48 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:46:48 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <string>
using namespace std;

namespace mu2e{
  namespace calorimeter{

class CrystalDetail{

public:

  CrystalDetail(
	       double crystalHalfTrans,
	       double crystalHalfLong,
	       string crystalMaterial,
	       string crystalWrapper,
	       double crystalWrapperThickness):
    _crystalHalfTrans       (crystalHalfTrans),
    _crystalHalfLong        (crystalHalfLong),
    _crystalMaterial        (crystalMaterial),
    _crystalWrapper         (crystalWrapper),
    _crystalWrapperThickness(crystalWrapperThickness) {}


  
  ~CrystalDetail () {};

  string crystalMaterial()         const { return _crystalMaterial;}
  string crystalWrapper()          const { return _crystalWrapper;}
  double crystalHalfTrans()        const { return _crystalHalfTrans; }
  double crystalHalfLong()         const { return _crystalHalfLong; }
  double crystalWrapperThickness() const { return _crystalWrapperThickness; }

private:

  double _crystalHalfTrans;
  double _crystalHalfLong;
  string _crystalMaterial;
  string _crystalWrapper;
  double _crystalWrapperThickness;


};


  } //namespace calorimeter
} //namespace mu2e

#endif
