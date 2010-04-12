
//
// Mu2e includes
#include "CalorimeterGeom/inc/Layer.hh"


using namespace std;

using CLHEP::Hep3Vector;


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
