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
#include "ConditionsService/inc/PhysicsParams.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {

  FoilParticleGenerator::FoilParticleGenerator(edm::RandomNumberGeneratorService::base_engine_t& engine,
                                               double tmin, double tmax):
    // time generation range
    _tmin ( tmin ),
    _tmax ( tmax ),
    // number of foils of the target
    _nfoils ( GeomHandle<Target>()->nFoils() ),
    // Random number distributions; getEngine comes from the base class.
    _randFlat ( engine ) ,
    _randTime ( engine ),
    _randFoils ( engine, &(binnedFoilsVolume()[0]), _nfoils )
  {
    // Check if nfoils is bigger than 0;
    if (_nfoils < 1) {
      throw cms::Exception("GEOM")
        << "no foils are present";
    }
  }
    
  FoilParticleGenerator::~FoilParticleGenerator()
  {
  }

  void FoilParticleGenerator::generatePositionAndTime(CLHEP::Hep3Vector& pos,
                                                      double& time) {

    // Get access to the geometry system.
    GeomHandle<Target> target;
    
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
    pos = CLHEP::Hep3Vector ( center.x()+r*cos(phi), 
                              center.y()+r*sin(phi), 
                              center.z()+dz );
    
    // Pick up a time random from am exponential distribution, with a given lifetime and in a defined range.
    ConditionsHandle<PhysicsParams> phyPar("ignored");
    double tau = phyPar->decayTime; 
    if (tau < 0 || tau > 3500) { //bigger than muon decay time
      throw cms::Exception("RANGE")
        << "nonsense decay time of bound state"; 
    }
    time = (_randTime.fire(_tmin, _tmax, tau));
  }
  
  
  vector<double> FoilParticleGenerator::binnedFoilsVolume() {
    
    vector<double> volumes;
    GeomHandle<Target> target;
    for (int i=0; i< _nfoils; ++i) {
      TargetFoil const& foil = target->foil(i);
      double rout = foil.rOut();
      double rin = foil.rIn();
      double halfthick = foil.halfThickness();
      double volume = CLHEP::pi*(rout-rin)*(rout-rin)*2.*halfthick;
      volumes.push_back(volume);
      // cout << "Foil " << i+1 << "  volume  " << volume << endl;
    }
    return volumes;
  } //FoilParticleGenerator::binnedFoilsVolume() 
  
}
