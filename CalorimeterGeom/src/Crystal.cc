//
// Hold information about one Straw.
//
// $Id: Crystal.cc,v 1.6 2010/05/17 21:47:33 genser Exp $
// $Author: genser $
// $Date: 2010/05/17 21:47:33 $
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
                CLHEP::Hep3Vector const& c,
                const CrystalDetail* detail,
                CLHEP::Hep3Vector const& w
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
