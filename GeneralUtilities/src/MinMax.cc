//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.cc,v 1.2 2011/12/16 23:12:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/12/16 23:12:30 $
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
