//
// enum to specify helicity, defined as the sign of the projection of a particles
// angular momentum onto its linear momentum vector
//  Original Author Dave Brown (LBNL) 15/8/2016
//
#ifndef DataProducts_Helicity_HH
#define DataProducts_Helicity_HH
#include <math.h>
#include <algorithm>
namespace mu2e {
  struct Helicity {
    enum helicity {neghel =-1, poshel =1, unknown=0};
    Helicity() { _value = unknown; }
    Helicity(int ival) { _value = static_cast<helicity>( (ival > 0) - (ival < 0) ); }
    Helicity(float val) { _value = val != 0.0 ? static_cast<helicity>(copysign(1.0,val)) : unknown; }
    bool operator == (Helicity const& other) const { return _value == other._value; }
    bool operator != (Helicity const& other) const { return !(operator ==(other)); }
    bool operator < (Helicity const& other) const { return _value < other._value; }
    helicity _value;
    static const char* name(Helicity const& h) {
      switch (h._value) {
	case Helicity::neghel:
	  return "Negative";
	  break;
	case Helicity::poshel:
	  return "Positive";
	  break;
	case Helicity::unknown : default:
	  return "Unknown";
	  break;
      }
    }
  };

}
#endif
