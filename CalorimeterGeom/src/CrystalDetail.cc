// $Id: CrystalDetail.cc,v 1.3 2010/04/13 17:23:48 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:23:48 $

// original authors Julie Managan and Robert Bernstein

//
// Mu2e includes
#include "CalorimeterGeom/inc/CrystalDetail.hh"

namespace mu2e{
  namespace calorimeter{

    CrystalDetail::CrystalDetail(
				 int    materialid,
				 double xhalfLength,
				 double yhalfLength,
				 double zhalfLength
				 ):
      _materialid(materialid),
      _xhalfLength(xhalfLength),
      _yhalfLength(yhalfLength),
      _zhalfLength(zhalfLength)
    {
    }
  
    CrystalDetail::~CrystalDetail (){
    }
  } //namespace calorimeter
} //namespace mu2e
