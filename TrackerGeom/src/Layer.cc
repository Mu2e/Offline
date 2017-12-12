//
// Hold information about one Layer in a tracker.
//
//
// $Id: Layer.cc,v 1.5 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Rob Kutschke
//

#include <sstream>
#include <iomanip>

#include "TrackerGeom/inc/Layer.hh"
#include "TrackerGeom/inc/Tracker.hh"

#ifndef __CINT__

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  Layer::Layer():
    _id(LayerId())
    // _nStraws(0)
    // _orig(CLHEP::Hep3Vector(0.,0.,0.)),
    // _delta(CLHEP::Hep3Vector(0.,0.,0.))
    {
  }

  // Layer::Layer(const LayerId& id
  //              // int      nStraws
  //              // const CLHEP::Hep3Vector& origin,
  //              // const CLHEP::Hep3Vector& delta
  //              ):
  //   _id(id)
  //   // _nStraws(nStraws)
  //   // _orig(origin),
  //   // _delta(delta)
  // {
  // }

  Layer::Layer(const LayerId& id ):
    _id(id){
  }

  string Layer::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getPlane() << "_"
       << _id.getPanel() << "_"
       << _id.getLayer();
    return os.str();
  }

  void Layer::fillPointers ( const Tracker& tracker ) const{

    // cout << __func__ << " Printing Layer nStraws before: " << nStraws() << endl;
    // for (int is = 0; is<(nStraws()*2); is+=2) {
    //   const Straw& straw = getStraw(is);
    //   cout << __func__ << " Printing Layer _straws info before: " << id()
    //        << setw(3) << is
    //        << " " << straw.id()
    //        << " " << straw.id2()
    //        << " " << straw.index()
    //        << endl;
    // }

    _straws.clear();
    for ( size_t i=0; i<_indices.size(); ++i ){
      StrawIndex idx = _indices[i];
      const Straw* straw =  &tracker.getStraw(idx);
      _straws.push_back(straw);
      straw->fillPointers(tracker);
    }
    _straw0Direction = _straws[0]->getDirection();
    _straw0MidPoint  = _straws[0]->getMidPoint();

    // cout << __func__ << " Printing Layer nStraws after: " << nStraws() << endl;
    // for (int is = 0; is<(nStraws()*2); is+=2) {
    //   const Straw& straw = getStraw(is);
    //   cout << __func__ << " Printing Layer _straws info after: " << id()
    //        << setw(3) << is
    //        << " " << straw.id()
    //        << " " << straw.id2()
    //        << " " << straw.index()
    //        << endl;
    // }

  }

} // namespace mu2e
#endif

