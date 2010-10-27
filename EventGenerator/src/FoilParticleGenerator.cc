//
// Generate a generic particle coming from target foils.
// Position, time and foil number of generated particle 
// is extracted random from appropiate distributions
//
// 
// Original author Gianni Onorato
// 
// 
// Notes
//
// 1) About the initialization of _randFoils.
//    The c'tor of RandGeneral wants, as its second argument, the starting
//    address of an array of doubles that describes the required shape.
//    The method binnedFoilsVolume returns, by value, a std::vector<double>.
//    We can get the required argument by taking the address of the first element 
//    of the std::vector. There is a subtlety about the return value of
//    those methods:  they return by value to a temporary variable that
//    we cannot see; this variable goes out of scope after the c'tor completes;
//    therefore its lifetime is managed properly.
//

// C++ includes.
#include <iostream>

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TMath.h"

using namespace std;


namespace mu2e {

  FoilParticleGenerator::FoilParticleGenerator(edm::RandomNumberGeneratorService::base_engine_t& engine ):
    // Random number distributions; getEngine comes from the base class.
    _randFlat ( engine ) ,
    _randTime ( engine ),
    _randFoils ( engine, &(binnedFoilsVolume()[0]), _nfoils )
  {
  }
    

  FoilParticleGenerator::~FoilParticleGenerator()
  {
  }

  void FoilParticleGenerator::generatePositionAndTime(CLHEP::Hep3Vector& pos,
                                                      double& time) {

    // Calculate spatial and temporal ranges
    _dcz  = (  _FPGczmax -  _FPGczmin);
    _dphi = ( _FPGphimax - _FPGphimin);
    _dt   = (   _FPGtmax -   _FPGtmin);
  
    
    // Get access to the geometry system.
    GeomHandle<Target> target;
    
    // Check if nfoils is bigger than 0;
    if (_nfoils < 1) {
      throw cms::Exception("GEOM")
        << "no foils are present";
    }
    
    // Pick a foil with a probability proportional to the volumes.
    int ifoil = static_cast<int>(_nfoils*_randFoils.fire());
    TargetFoil const& foil = target->foil(ifoil);
    
    // Foil properties.
    CLHEP::Hep3Vector const& center = foil.center();
    const double r1 = foil.rIn();
    const double dr = foil.rOut() - r1;
    
    // A random point within the foil. To change.
    const double r   = r1 + dr*_randFlat.fire();
    const double dz  = (-1.+2.*_randFlat.fire())*foil.halfThickness();
    const double phi = CLHEP::twopi*_randFlat.fire();
    CLHEP::Hep3Vector tpos( center.x()+r*cos(phi), 
                             center.y()+r*sin(phi), 
                             center.z()+dz );

    pos = tpos;

    // Pick up a time random from am exponential distribution. To improve
    double ttime = _FPGtmin   +  _dt*(_randTime.fire(1/(_FPGtmax)));
    time = ttime;
  }
  
  
  vector<double> FoilParticleGenerator::binnedFoilsVolume() {
    
    vector<double> volumes;
    GeomHandle<Target> target;
    _nfoils = target->nFoils();
    for (int i=0; i< _nfoils; i++) {
      TargetFoil const& foil = target->foil(i);
      double rout = foil.rOut();
      double rin = foil.rIn();
      double halfthick = foil.halfThickness();
      double volume = CLHEP::pi*(rout-rin)*(rout-rin)*2*halfthick;
      volumes.push_back(volume);
      // cout << "Foil " << i+1 << "  volume  " << volume << endl;
    }
    return volumes;
  } //FoilParticleGenerator::binnedFoilsVolume() 
  
}
