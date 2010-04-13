#ifndef CRYSTALDETAIL_HH
#define CRYSTALDETAIL_HH
// $Id: CrystalDetail.hh,v 1.2 2010/04/13 17:14:51 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:14:51 $

// original authors Julie Managan and Robert Bernstein

namespace mu2e{
  namespace calorimeter{

class CrystalDetail{

public:
  CrystalDetail():
    _materialid(-1),
    _xhalfLength(-1),
    _yhalfLength(-1),
    _zhalfLength(-1)
  {}

  CrystalDetail( int materialid,
	       double xhalfLength,
	       double yhalfLength,
	       double zhalfLength
	       );

  
  ~CrystalDetail ();

  int materialId() const { return _materialid;}

  double xhalfLength() const { return _xhalfLength; }
  double yhalfLength() const { return _yhalfLength; }
  double zhalfLength() const { return _zhalfLength; }

  // Compiler generated copy and assignment constructors
  // should be OK.
  
private:
  int    _materialid;
  double _xhalfLength;
  double _yhalfLength;
  double _zhalfLength;

};


  } //namespace calorimeter
} //namespace mu2e

#endif
