#ifndef GeneralUtilities_MinMax_hh
#define GeneralUtilities_MinMax_hh

//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.hh,v 1.6 2011/12/16 23:12:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/12/16 23:12:30 $
//
// Original author Rob Kutschke

#include <limits>
#include <iostream>

class MinMax{

public:

  MinMax():
    _n(0),
    _min(   std::numeric_limits<double>::max()),
    _max(  -std::numeric_limits<double>::max()),
    _small( std::numeric_limits<double>::max())
  {}

  MinMax(double x):
    _n(1),
    _min(   std::numeric_limits<double>::max()),
    _max(  -std::numeric_limits<double>::max()),
    _small( std::numeric_limits<double>::max()){
    accumulate(x);
  }

  // Accept compiler written d'tor, copy c'tor and assignment operator.

  // Update min/max/smallest with this entry.
  void accumulate(double x);

  // Accessors
  double n()        const { return _n;}
  double min()      const { return _min;}
  double max()      const { return _max;}
  double smallest() const { return _small;}
  double delta()    const { return _max-_min; }

private:

  // Limiting values of the numbers compared so far.
  int    _n;
  double _min;
  double _max;
  double _small;

};

inline std::ostream& operator<<(std::ostream& ost,
                               const MinMax&  mm ){
  ost << "( "
      << mm.n()        << " "
      << mm.min()      << " "
      << mm.max()      << " "
      << mm.smallest()
      << " )";
  return ost;
}

#endif /* GeneralUtilities_MinMax_hh */
