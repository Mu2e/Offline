//
// Hold information about one Straw.
//
// $Id: Crystal.cc,v 1.5 2010/05/12 14:58:49 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/05/12 14:58:49 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <sstream>
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#ifndef __CINT__ 

using std::vector;
using CLHEP::Hep3Vector;

namespace mu2e {
  namespace calorimeter{
  Crystal::Crystal( const CrystalId& id,
                CrystalIndex index,
                Hep3Vector const& c,
                const CrystalDetail* detail,
                Hep3Vector const& w
                ):
    _id(id),
    _index(index),
    _c(c),
    _detail(detail),
    _w(w.unit()){
  }
  } //namespace calorimeter
} // namespace mu2e
#endif
