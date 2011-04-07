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
#include <fstream>

//Framework includes
#include "FWCore/ParameterSet/interface/FileInPath.h"

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

static const double timeMaxDelay = 3000;
static const int nBinsForTimeDelayPDF = 150;

namespace mu2e {
  
  FoilParticleGenerator::FoilParticleGenerator(edm::RandomNumberGeneratorService::base_engine_t& engine,
                                               double tmin, double tmax, 
                                               foilGen_enum foilAlgo, 
                                               posGen_enum  posAlgo, 
                                               timeGen_enum  timeAlgo,
                                               bool PTtoSTdelay):
    // time generation range
    _tmin ( tmin ),
    _tmax ( tmax ),
    // selected algorithm for foils 
    _foilAlgo ( foilAlgo ),
    _posAlgo ( posAlgo ),
    _timeAlgo ( timeAlgo ),
    // number of foils of the target
    _nfoils ( GeomHandle<Target>()->nFoils() ),
    // Random number distributions; getEngine comes from the base class.
    _randFlat ( engine ) ,
    _randTime( engine ),
    _randFoils ( engine, &(binnedFoilsVolume()[0]), _nfoils ),
    _randExpoFoils ( engine, &(weightedBinnedFoilsVolume()[0]), _nfoils ),
    _delayTime( engine, &(timePathDelay()[0]), nBinsForTimeDelayPDF ), 
    _PTtoSTdelay ( PTtoSTdelay )
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

    
    // Pick a foil 
    
    switch (_foilAlgo) {
    case flatFoil:
      _ifoil = getFlatRndFoil();
      break;
    case volWeightFoil:
      _ifoil = getVolumeRndFoil();
      break;
    case expoVolWeightFoil:
      _ifoil = getVolumeAndExpoRndFoil();
      break;
    default:
      break;
    }

    // Get access to the geometry system.
    GeomHandle<Target> target;
    TargetFoil const& foil = target->foil(_ifoil);

    //Pick up position
    switch (_posAlgo) {
    case flatPos:
      pos = getFlatRndPos(foil);
      break;
    default:
      break;
    }

    //Pick up time
    switch (_timeAlgo) {
    case flatTime:
      time = getFlatRndTime();
      break;
    case limitedExpoTime:
      time = getLimitedExpRndTime();
      break;
    default:
      break;
    }

    if (_PTtoSTdelay) {
      double deltat = includeTimeDelay();
      time += deltat;
    }
    
  }
  

  int FoilParticleGenerator::iFoil() {
    return _ifoil;
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
  
  
  //For Pi Capture production: the previous code used a randexponential to describe
  // the generation in foils. Lambda of the distribution was 1. Since the dist
  // output goes over 1, it was forced to regenerate the rnd number if bigger than 1.
  //Now I include thisproduction in the "foil volume weighted" frame.
  //The 0-1 output of the exponential generation is divided in (_nfoils) bins.
  //The integral of the exponential between i/_nfoils and (i+1)/_nfoils, divided for
  //the x-axis step of the integral (1/_nfoils), is the mean value of the exponential
  //function in the bin corresponding to a foil. I use this value as a weight for
  //the foil volume associated to the bin. 
  //Procedure surely to refine.
  
  vector<double> FoilParticleGenerator::weightedBinnedFoilsVolume() {
    
    vector<double> volumes = binnedFoilsVolume();
    if (volumes.size()!= (size_t) _nfoils) {
      throw cms::Exception("GEOM")
        << "something wrong in number of foils";
    }
    double step = 1./_nfoils;
    for (int i=0; i< _nfoils; ++i) {
      cout << volumes[i] << '\t';
      double weight = (exp(-(step*i))*(1-(exp(-step))))/step;
      cout << weight << '\t';
      volumes[i] = volumes[i]*weight;
      cout << volumes[i] << endl;
    }
    return volumes;
  } //FoilParticleGenerator::weightedBinnedFoilsVolume() 
  

  vector<double> FoilParticleGenerator::timePathDelay() {

    vector<double> muonTimeDelay;
    edm::FileInPath muonDelayFileName("ConditionsService/data/timeDelayDist.txt");
    string MuonFileFIP = muonDelayFileName.fullPath();
    fstream infile(MuonFileFIP.c_str(), ios::in);
    if (infile.is_open()) {
      double val;
      
      for (int i=0; i < nBinsForTimeDelayPDF; ++i) {
        infile >> val;
        muonTimeDelay.push_back(val);
      }
    } else {
      cout << "No file associated for the muon arriving delay distribution" << endl;
      for (int i=0; i < nBinsForTimeDelayPDF; ++i) {
        muonTimeDelay.push_back(1);
      }
    }
    
    return muonTimeDelay;
    
  }


  // Pick up a random foil from a flat distribution
  int FoilParticleGenerator::getFlatRndFoil() {
    return static_cast<int>(_nfoils*_randFlat.fire());
  }

  // Pick up a random foil from a flat distribution 
  // weighted by foil volume
  int FoilParticleGenerator::getVolumeRndFoil() {
    return  static_cast<int>(_nfoils*_randFoils.fire());
  }

  // Pick up a random foil from a negative exponential 
  // distribution weighted by foil volume 
  int FoilParticleGenerator::getVolumeAndExpoRndFoil() {
    return static_cast<int>(_nfoils*_randExpoFoils.fire());
  }

  // Pick up a random position within the foil
  CLHEP::Hep3Vector FoilParticleGenerator::getFlatRndPos(TargetFoil const& theFoil) {
    
    // Foil properties.
    CLHEP::Hep3Vector const& center = theFoil.center();
    const double r1 = theFoil.rIn();
    const double dr = theFoil.rOut() - r1;
    
    // A random point within the foil.
    const double r   = r1 + dr*_randFlat.fire();
    const double dz  = (-1.+2.*_randFlat.fire())*theFoil.halfThickness();
    const double phi = CLHEP::twopi*_randFlat.fire();
    return CLHEP::Hep3Vector ( center.x()+r*cos(phi), 
                               center.y()+r*sin(phi), 
                               center.z()+dz );

  }

  // Pick up a random generation time from a flat distribution
  double FoilParticleGenerator::getFlatRndTime() {
    return _tmin + _randFlat.fire() * (_tmax-_tmin);
  }

  double FoilParticleGenerator::includeTimeDelay() {

    double dt = timeMaxDelay * _delayTime.fire();
    return dt;
  }


  // Pick up a time random from am exponential distribution, 
  // with a given lifetime and in a defined range.
  double FoilParticleGenerator::getLimitedExpRndTime() {
    ConditionsHandle<PhysicsParams> phyPar("ignored");
    double tau = phyPar->decayTime; 
    if (tau < 0 || tau > 3500) { //bigger than muon decay time
      throw cms::Exception("RANGE")
        << "nonsense decay time of bound state"; 
    }
    return _randTime.fire(_tmin, _tmax, tau);
    
  }
}


