//
// G4 interface to the Detector Solenoid full magnetic field.
//
// $Id: DSField.cc,v 1.3 2010/05/18 22:33:49 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:33:49 $
//
// Original author Julie Managan and Bob Bernstein

#include <cmath>
#include <iostream>

#include "Mu2eG4/inc/DSField.hh"
#include "Mu2eG4/inc/BFwcont.hh"
#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using namespace std;

namespace mu2e {

  void DSField::GetFieldValue(const G4double Point[4],
                              G4double *Bfield) const {

    CLHEP::Hep3Vector point(Point[0],Point[1],Point[2]);
    CLHEP::Hep3Vector BField = _p->GetBField(point);
    Bfield[0] = BField.x()*tesla; 
    Bfield[1] = BField.y()*tesla; 
    Bfield[2] = BField.z()*tesla;

    //  cout << "(" << BField.x() << "," << BField.y() << "," << BField.z() << ")" << endl; 

  }

} // end namespace mu2e
