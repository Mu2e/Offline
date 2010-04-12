
//
// Mu2e includes
#include "Calorimeter/inc/CrystalDetail.hh"


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
