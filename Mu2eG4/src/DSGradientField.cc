//
// G4 interface to the Detector Solenoid gradient magnetic field.
//
// Author Ivan Logashenko
//
// C++ includes
//#include <cmath>
#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/DSGradientField.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using namespace std;

namespace mu2e {

  DSGradientField::DSGradientField( std::string name, G4ThreeVector mapOrigin,
                                    double gradient, double B0 ):
    _name(name), _mapOrigin(mapOrigin), _gradient(gradient), _B0(B0) {
  }

  // This is the entry point called by G4.
  void DSGradientField::GetFieldValue(const G4double Point[4],
                                      G4double *Bfield) const {
      
//      cout << "calling DSGradientField::GetFieldValue()" << endl;


    CLHEP::Hep3Vector point(Point[0],Point[1],Point[2]);

    Bfield[2] = _B0 + _gradient*(point-_mapOrigin).z();

    //double r = (point-_mapOrigin).perp();
    //double Br = - _gradient*r/2.0;

    //Bfield[0] = (point.x()-_mapOrigin.x())/r*Br;
    Bfield[0] = -(point.x()-_mapOrigin.x())*_gradient/2.0;
    Bfield[1] = -(point.y()-_mapOrigin.y())*_gradient/2.0;
      
//      cout << "Mu2eGlobalField map, "
//      << " point=("<<point.x()<<","<<point.y()<<","<<point.z()<<")"
//      << " field=("<<Bfield[0]<<","<<Bfield[1]<<","<<Bfield[2]<<")"
//      << endl;


  }

} // end namespace mu2e
