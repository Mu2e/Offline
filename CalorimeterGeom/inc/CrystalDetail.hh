#ifndef CRYSTALDETAIL_HH
#define CRYSTALDETAIL_HH
// $Id: CrystalDetail.hh,v 1.6 2010/05/18 22:07:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:07:16 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <string>

namespace mu2e{
  namespace calorimeter{

class CrystalDetail{

public:

  CrystalDetail(
               double      crystalHalfTrans,
               double      crystalHalfLong,
               std::string crystalMaterial,
               std::string crystalWrapper,
               double      crystalWrapperHalfThickness):
    _crystalHalfTrans            (crystalHalfTrans),
    _crystalHalfLong             (crystalHalfLong),
    _crystalMaterial             (crystalMaterial),
    _crystalWrapper              (crystalWrapper),
    _crystalWrapperHalfThickness (crystalWrapperHalfThickness) {}


  
  ~CrystalDetail () {};

  const std::string& getCrystalMaterial()         const{ return _crystalMaterial;}
  const std::string& getCrystalWrapper()          const{ return _crystalWrapper;}
  double  getCrystalHalfTrans()              const{ return _crystalHalfTrans; }
  double  getCrystalHalfLong()               const{ return _crystalHalfLong; }
  double  getCrystalWrapperHalfThickness()   const{ return _crystalWrapperHalfThickness; }

private:

  double      _crystalHalfTrans;
  double      _crystalHalfLong;
  std::string _crystalMaterial;
  std::string _crystalWrapper;
  double      _crystalWrapperHalfThickness;


};


  } //namespace calorimeter
} //namespace mu2e

#endif
