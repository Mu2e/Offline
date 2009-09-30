#ifndef STRAWDETAIL_HH
#define STRAWDETAIL_HH

//
// Incomplete prototype of a class to hold information about
// the properties of each type of straw.  We will 
// definitely need different types of straws for the 
// conducting and non-conducting straws.
// 

//
// $Id: StrawDetail.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//


#include <vector>
#include <string>

namespace mu2e {

class StrawDetail{

public:
  StrawDetail():
    _materialNames(),
    _radius(-1.),
    _thickness(-1.),
    _halfLength(-1.),
    _rwire(-1.)
  {}

  StrawDetail( std::vector<std::string> const& materialNames,
	       double radius,
	       double thickness,
	       double length,
	       double rwire
	       );
  
  ~StrawDetail ();

  std::string const& materialName(int idx) const { return _materialNames.at(idx);}

  std::vector<std::string> materialNames() const { return _materialNames;}

  double      outerRadius()   const { return _radius;}
  double      innerRadius()   const { return _radius-_thickness;}
  double      thickness()     const { return _thickness; }
  double      length()        const { return _halfLength*2.0; }
  double      halfLength()    const { return _halfLength; }

  // Compiler generated copy and assignment constructors
  // should be OK.
  
private:

  // Order of materials is:
  // Straw itself, Gas volume, wire;
  std::vector<std::string> _materialNames;

  double _radius;
  double _thickness;
  double _halfLength;
  double _rwire;

};

}  //namespace mu2e

#endif
