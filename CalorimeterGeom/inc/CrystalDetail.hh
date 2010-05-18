#ifndef CRYSTALDETAIL_HH
#define CRYSTALDETAIL_HH
// $Id: CrystalDetail.hh,v 1.5 2010/05/18 20:29:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 20:29:09 $

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
               double crystalWrapperHalfThickness):
    _crystalHalfTrans            (crystalHalfTrans),
    _crystalHalfLong             (crystalHalfLong),
    _crystalMaterial             (crystalMaterial),
    _crystalWrapper              (crystalWrapper),
    _crystalWrapperHalfThickness (crystalWrapperHalfThickness) {}


  
  ~CrystalDetail () {};

  const string& getCrystalMaterial()         const{ return _crystalMaterial;}
  const string& getCrystalWrapper()          const{ return _crystalWrapper;}
  double  getCrystalHalfTrans()              const{ return _crystalHalfTrans; }
  double  getCrystalHalfLong()               const{ return _crystalHalfLong; }
  double  getCrystalWrapperHalfThickness()   const{ return _crystalWrapperHalfThickness; }

private:

  double _crystalHalfTrans;
  double _crystalHalfLong;
  string _crystalMaterial;
  string _crystalWrapper;
  double _crystalWrapperHalfThickness;


};


  } //namespace calorimeter
} //namespace mu2e

#endif
