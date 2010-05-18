#include <cmath>
#include <iostream>

#include "Mu2eG4/inc/DSField.hh"
//#include "../inc/BFBMvec.hh"
#include "Mu2eG4/inc/BFwcont.hh"
#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using namespace std;

void DSField::GetFieldValue(const G4double Point[4],
                             G4double *Bfield) const {

  CLHEP::Hep3Vector point(Point[0],Point[1],Point[2]);
  CLHEP::Hep3Vector BField = _p->GetBField(point);
  Bfield[0] = BField.x()*tesla; 
  Bfield[1] = BField.y()*tesla; 
  Bfield[2] = BField.z()*tesla;

  //  cout << "(" << BField.x() << "," << BField.y() << "," << BField.z() << ")" << endl; 

}

