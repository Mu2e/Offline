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
#include "TrkExt/inc/TrkExtDetectors.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtDetectors::TrkExtDetectors() : TrkExtShape()
  { 
  }

  void TrkExtDetectors::initialize() {
    _ds.initialize();
    _pa.initialize();
    _st.initialize();
  }

  bool TrkExtDetectors::contains (Hep3Vector &xx) {
    return _ds.contains(xx);
  }

  TrkExtDetectorList::Enum TrkExtDetectors::volumeId(Hep3Vector & xx) {
    if      (_pa.contains(xx)) return TrkExtDetectorList::ProtonAbsorber;
    else if (_st.contains(xx)) return TrkExtDetectorList::StoppingTarget;
    else if (_ds.contains(xx)) return TrkExtDetectorList::ToyDS;
    else                       return TrkExtDetectorList::Undefined;
  }


    
  Hep3Vector  TrkExtDetectors::intersection (const Hep3Vector & x1, const Hep3Vector & x2) {
    TrkExtDetectorList::Enum f1, f2;
    Hep3Vector xstart = x1;
    Hep3Vector xstop = x2;

    f1 = volumeId(xstart);
    f2 = volumeId(xstop);
    if (f1 == f2) {
      cerr << "TrkExtDetectors Warning : call intersection at wrong positions at begin " << f1 << xstart << ", " << f2 <<xstop<< endl;
      Hep3Vector ret = (x1+x2)*0.5;
      return ret;
    }
    int i = 0;
    
    do {

      Hep3Vector xhalf = (xstart+xstop)*0.5;
      TrkExtDetectorList::Enum fhalf = volumeId(xhalf);
      if (fhalf == f1) {
        xstart = xhalf;
      }
      else if (fhalf == f2) {
        xstop = xhalf;
      }
      else {
        cerr << "TrkExtDetectors Warning : call intersection at wrong positions at processing "<< f1 << ", " << f2  << endl;
        return xhalf;
      }
      ++i;
    } while ((xstart-xstop).mag() > _limit);
    return xstop;
  }

  double TrkExtDetectors::mostProbableEnergyLoss (const Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid) {
    switch (volid) {
      case TrkExtDetectorList::ProtonAbsorber:
        return _pa.mostProbableEnergyLoss(p, ds);
      case TrkExtDetectorList::StoppingTarget:
        return _st.mostProbableEnergyLoss(p, ds);
      case TrkExtDetectorList::Undefined:
      case TrkExtDetectorList::ToyDS:
      default:
        return 0;
    }
    return 0;
  }


  double TrkExtDetectors::meanEnergyLoss (const Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid) {
    switch (volid) {
      case TrkExtDetectorList::ProtonAbsorber:
        return _pa.meanEnergyLoss(p, ds);
      case TrkExtDetectorList::StoppingTarget:
        return _st.meanEnergyLoss(p, ds);
      case TrkExtDetectorList::Undefined:
      case TrkExtDetectorList::ToyDS:
      default:
        return 0;
    }
    return 0;
  }

  double TrkExtDetectors::scatteringAngle (const CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid) {
    switch (volid) {
      case TrkExtDetectorList::ProtonAbsorber:
        return _pa.scatteringAngle(p, ds);
      case TrkExtDetectorList::StoppingTarget:
        return _st.scatteringAngle(p, ds);
      case TrkExtDetectorList::Undefined:
      case TrkExtDetectorList::ToyDS:
      default:
        return 0;
    }
    return 0;
  }



} // end namespace mu2e

