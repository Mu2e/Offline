
//
// C++ includes
#include <cmath>
#include <cstdlib>

//
// Mu2e includes
#include "Calorimeter/inc/Crystal.hh"


//
// other includes
#include "CLHEP/Vector/ThreeVector.h"
#ifndef __CINT__ 

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

#endif
