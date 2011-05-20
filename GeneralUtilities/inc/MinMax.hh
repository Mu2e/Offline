#ifndef GeneralUtilities_MinMax_hh
#define GeneralUtilities_MinMax_hh

//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.hh,v 1.5 2011/05/20 20:18:23 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:18:23 $
//
// Original author Rob Kutschke

#include <cmath>
#include <limits>

class MinMax{

public:

  MinMax():
    _min(   std::numeric_limits<double>::max()),
    _max(  -std::numeric_limits<double>::max()),
    _small( std::numeric_limits<double>::max())
  {}

  MinMax(double x):
    _min(   std::numeric_limits<double>::max()),
    _max(  -std::numeric_limits<double>::max()),
    _small( std::numeric_limits<double>::max()){
    compare(x);
  }

  // Compare this x to the previously existing min/max
  void compare(double x);

  // Return limiting values.
  double min() const {return _min;}
  double max() const {return _max;}
  double smallest() const {return _small;}

private:

  // Limiting values of the numbers compared so far.
  double _min;
  double _max;
  double _small;

};

#endif /* GeneralUtilities_MinMax_hh */
