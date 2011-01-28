#ifndef MinMax_HH
#define MinMax_HH

//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.hh,v 1.2 2011/01/28 23:51:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/01/28 23:51:57 $
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

#endif
