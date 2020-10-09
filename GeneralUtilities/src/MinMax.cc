//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
//
// Original author Rob Kutschke

#include <limits>
#include <cmath>

#include "GeneralUtilities/inc/MinMax.hh"

using namespace std;

void MinMax::accumulate(double x){
  ++_n;
  _min   = (x < _min ) ? x : _min;
  _max   = (x > _max ) ? x : _max;
  _small = ( abs(x) < _small ) ? abs(x) : _small;
}
