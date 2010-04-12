#ifndef CRYSTALDETAIL_HH
#define CRYSTALDETAIL_HH

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

#endif
