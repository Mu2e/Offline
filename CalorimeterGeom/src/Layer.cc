// $Id: Layer.cc,v 1.3 2010/04/13 17:23:33 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:23:33 $

// original authors Julie Managan and Robert Bernstein

//
// Mu2e includes
#include "CalorimeterGeom/inc/Layer.hh"


using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e{
  namespace calorimeter{

Layer::Layer():
  _id(LayerId()),
  _nCrystals(0),
  _orig(Hep3Vector(0.,0.,0.)),
  _delta(Hep3Vector(0.,0.,0.))
{
}

Layer::Layer(const LayerId& id,
                   int      nCrystals,
	     const Hep3Vector& origin,
	     const Hep3Vector& delta
	     ):
  _id(id),
  _nCrystals(nCrystals),
  _orig(origin),
  _delta(delta)
{
}

Layer::Layer(const LayerId& id ):
  _id(id)
{
}

Layer::~Layer(){}

  } //namespace calorimeter
} //namespace mu2e
