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

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

static mu2e::ConfigFileLookupPolicy findConfig;
static const string StMuFileString = findConfig("ConditionsService/data/StoppedMuons.txt");
static const double timeMaxDelay = 3000;
static const double nBinsForTimeDelayPDF = 150;
static fstream inMuFile(StMuFileString.c_str(), ios::in);

namespace mu2e {

  FoilParticleGenerator::FoilParticleGenerator(art::RandomNumberGenerator::base_engine_t& engine,
                                               double tmin, double tmax,
                                               foilGen_enum foilAlgo,
                                               posGen_enum  posAlgo,
                                               timeGen_enum  timeAlgo,
                                               bool targetFrame,
                                               bool PTtoSTdelay,
                                               bool pPulseDelay ):
    _prodTargetOffset(),
    _prodTargetCenter(),
    _g4beamlineOrigin(),
    _g4beamlineExtraOffset(CLHEP::Hep3Vector()),
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
    _muTimeDecay (getMuTimeDecay()),
    _randNegExpoTime( engine, -_muTimeDecay ),
    _randFoils ( engine, &(binnedFoilsVolume()[0]), _nfoils ),
    _randExpoFoils ( engine, &(weightedBinnedFoilsVolume()[0]), _nfoils ),
    _delayTime( engine, &(timePathDelay()[0]), nBinsForTimeDelayPDF ),
    _pulseTime( engine ),
    _targetFrame( targetFrame ),
    _PTtoSTdelay ( PTtoSTdelay ),
    _pPulseDelay ( pPulseDelay )
  {


    CLHEP::Hep3Vector offset_default(0.0,0.0,1764.5);
    _prodTargetOffset = offset_default;

    CLHEP::Hep3Vector g4beamlineOrigin_default(3904.,0.,-7929.);
    _g4beamlineOrigin = g4beamlineOrigin_default;

    CLHEP::Hep3Vector g4beamlineExtraOffset_default;
    _g4beamlineExtraOffset = g4beamlineExtraOffset_default;

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const& geomConfig = geom->config();
    _prodTargetCenter = geomConfig.getHep3Vector("productionTarget.position");;

  // Check if nfoils is bigger than 0;
    if (_nfoils < 1) {
      throw cet::exception("GEOM")
        << "no foils are present";
    }

    if (!inMuFile.is_open()) {
      throw cet::exception("GEOM")
        << "no stopped muon file is present";
    }

  }

  FoilParticleGenerator::~FoilParticleGenerator()
  {
  }

  void FoilParticleGenerator::generatePositionAndTime(CLHEP::Hep3Vector& pos,
                                                      double& time) {


    // Pick a foil

    time = -1000;
    while (time < _tmin) {

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
      case muonFileInputFoil:
        _ifoil = 0;
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
      case muonFileInputPos:
        getInfoFromFile(pos, time);
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
      case negExp:
        time += getNegativeExpoRndTime();
      default:
        break;
      }

      if (_pPulseDelay) {
	pulseDelay = includePulseDelay();
	time += pulseDelay;
      } else {
	pulseDelay = 0.;
      }

      if (_PTtoSTdelay) {
        muDelay = includeTimeDelay();
        time += muDelay;
      } else {
	muDelay = 0.;
      }



    }
    if (_posAlgo==muonFileInputPos) {
      if( _targetFrame ) {
        pos -= _prodTargetOffset;           // Move to target coordinate system
        pos += _prodTargetCenter; // Move to Mu2e coordinate system
      } else{
        pos += _g4beamlineOrigin;
        pos += _g4beamlineExtraOffset;
      }
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
      throw cet::exception("GEOM")
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
    string MuonFileFIP = findConfig("ConditionsService/data/timeDelayDist.txt");
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

  //Pick up the position from the input stopped muon file
  void FoilParticleGenerator::getInfoFromFile(CLHEP::Hep3Vector& pos, double& time) {

    //Start reading the input file from the beginning when it reaches the end
    bool gotthem = false;
    string line;
    while (!gotthem) {
      while(getline(inMuFile,line)) {
        stringstream s_line;
        s_line << line;
        size_t skip = line.find("\\\\");
        if (skip != string::npos) continue;
        if (line == "") continue;
        double x, y, z, t;
        s_line >> x >> y >> z >> t;
        CLHEP::Hep3Vector temppos(x,y,z);
        pos = temppos;
        time = t;
        gotthem = true;
        return;
      }
      inMuFile.clear();
      inMuFile.seekg(0, ios::beg);
    }

    return;

  }

  //Add a time interval taken from a negative exponential pdf
  double FoilParticleGenerator::getNegativeExpoRndTime() {
    return _randNegExpoTime.fire();
  }


  double FoilParticleGenerator::getMuTimeDecay() {

  ConditionsHandle<PhysicsParams> phyPar("ignored");
  double tau = phyPar->decayTime;
  if (tau < 0 || tau > 3500) { //bigger than muon decay time
    throw cet::exception("RANGE")
      << "nonsense decay time of bound state";
    }
  return tau;
  }

  // Pick up a random generation time from a flat distribution
  double FoilParticleGenerator::getFlatRndTime() {
    return _tmin + _randFlat.fire() * (_tmax-_tmin);
  }

  double FoilParticleGenerator::includeTimeDelay() {

    double dt = timeMaxDelay * _delayTime.fire();
    return dt;
  }


  double FoilParticleGenerator::includePulseDelay() {

  double dt = _pulseTime.fire();
  return dt;

  }


  // Pick up a time random from am exponential distribution,
  // with a given lifetime and in a defined range.
  double FoilParticleGenerator::getLimitedExpRndTime() {

    return _randTime.fire(0, _tmax, _muTimeDecay);

  }
}


