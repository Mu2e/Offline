#ifndef GeneralUtilities_MinMax_hh
#define GeneralUtilities_MinMax_hh

//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.hh,v 1.3 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke

#include <limits>
#include <cmath>

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
