//
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtShape.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtShape::TrkExtShape( double boundaryLimit ) :
    _limit(boundaryLimit)
  {  }

  Hep3Vector  TrkExtShape::intersection (const Hep3Vector & x1, const Hep3Vector & x2) {
    cout << "TrkExtShape called" << endl;
    bool f1, f2;
    Hep3Vector xstart = x1;
    Hep3Vector xstop = x2;
    do {
      f1 = contains(xstart);
      f2 = contains(xstop);
      if (f1 == f2) {
        cerr << "TrkExtShape Warning : call intersection at wrong positions at begin" << endl;
        Hep3Vector ret = (x1+x2)*0.5;
        return ret;
      }

      Hep3Vector xhalf = (xstart+xstop)*0.5;
      bool fhalf = contains(xhalf);
      if (fhalf == f1) {
        xstop = xhalf;
      }
      else if (fhalf == f2) {
        xstart = xhalf;
      }
      else {
        cerr << "TrkExtShape Warning : call intersection at wrong positions at processing" << endl;
        return xhalf;
      }
    } while ((xstart-xstop).mag() > _limit);
    return xstop;
  }





} // end namespace mu2e

