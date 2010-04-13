// $Id: Crystal.cc,v 1.3 2010/04/13 17:23:41 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:23:41 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <cmath>
#include <cstdlib>

//
// Mu2e includes
#include "CalorimeterGeom/inc/Crystal.hh"


//
// other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e{
  namespace calorimeter{

    using CLHEP::Hep3Vector;

    Crystal::Crystal():
      _id(CrystalId()),
      _index(0),
      _cNominal(Hep3Vector(0.,0.,0.)),
      _detail(0),
      _detailIndex(0),
      _w(Hep3Vector(0.,0.,1.))
    {
    }


    Crystal::Crystal( const CrystalId& id,
		      CrystalIndex  index,
		      Hep3Vector const& c,
		      const CrystalDetail* detail,
		      int detailIndex,
		      double wtx,
		      double wty
		      ):
      _id(id),
      _index(index),
      _cNominal(c),
      _detail(detail),
      _detailIndex(detailIndex)
    {
      _w = Hep3Vector(wtx,wty,1.).unit();
    }

    Crystal::Crystal( const CrystalId& id,
		      CrystalIndex index,
		      Hep3Vector const& c,
		      const CrystalDetail* detail,
		      int detailIndex,
		      Hep3Vector const& w
		      ):
      _id(id),
      _index(index),
      _cNominal(c),
      _detail(detail),
      _detailIndex(detailIndex),
      _w(w.unit())
    {
    }

    Crystal::~Crystal (){
    }

  } //namespace calorimeter
} //namespace mu2e

