//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.cc,v 1.1 2010/06/22 16:05:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:05:18 $
//
// Original author Rob Kutschke

#include <limits>
#include <cmath>

#include "GeneralUtilities/inc/MinMax.hh"

using namespace std;

void MinMax::compare(double x){
  _min   = (x < _min ) ? x : _min;
  _max   = (x > _max ) ? x : _max;
  _small = ( abs(x) < _small ) ? abs(x) : _small;
}
