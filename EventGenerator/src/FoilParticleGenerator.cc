//
// Generate a generic particle coming from target foils.
// Position, time and foil number of generated particle 
// is extracted random from appropiate distributions
//
// 
// Original author Gianni Onorato
// 
// Notes
// 1) This code uses a incorrect model of the distribution of DIO's over the
//    targets.  It is uniform across each target.
//    At a future date this needs to be made more realistic.
// 2) This code uses an incorrect model of the distribution of DIO's in time.
//    At a future date this needs to be made more realistic.
// 3) This codes uses (Emax-E)**5 for the momentum distribution.  At a future
//    date this needs to be improved.
// 4) About the initialization of _shape and _randFoils.
//    The c'tor of RandGeneral wants, as its second argument, the starting
//    address of an array of doubles that describes the required shape.
//    The methods binnedEnergySpectrum and binnedFoilsVolume return, by value, a std::vector<double>.
//    We can get the required argument by taking the address of the first element 
//    of the std::vector. There is a subtlety about the return value of
//    those methods:  they return by value to a temporary variable that
//    we cannot see; this variable goes out of scope after the c'tor completes;
//    therefore its lifetime is managed properly.
//

// C++ includes.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

// General Utilities
#include "GeneralUtilities/inc/pow.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;


namespace mu2e {

  FoilParticleGenerator::FoilParticleGenerator()
  {
  }

  FoilParticleGenerator::FoilParticleGenerator( PDGCode::type pdgId, 
                                                GenId::enum_type genId):
    _pdgId(pdgId),
    _genId(genId) {

    // Get the particle mass from the particle data table (in MeV).
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    const HepPDT::ParticleData& part_data = pdt->particle(_pdgId);
    _mass = part_data.mass().value();
    _nfoils = 0;
  }
  
  FoilParticleGenerator::~FoilParticleGenerator() 
  {
  }
  
  void FoilParticleGenerator::setRandomEngine(edm::RandomNumberGeneratorService::base_engine_t& engine ) {

    // Random number distributions; getEngine comes from the base class.
    _randFlat = auto_ptr<CLHEP::RandFlat> ( new CLHEP::RandFlat( engine ) ); 
    _randFoils = auto_ptr<CLHEP::RandGeneral> ( new CLHEP::RandGeneral( engine, &(binnedFoilsVolume()[0]), _nfoils) ); 
    _randomUnitSphere = auto_ptr<RandomUnitSphere> ( new RandomUnitSphere( engine, _czmin, _czmax, _phimin, _phimax) );
    _randTime = auto_ptr<CLHEP::RandExponential> (new CLHEP::RandExponential( engine ) );
    if (_genId == GenId::dio1) {
      // See Note 4.
      _shape = auto_ptr<CLHEP::RandGeneral> ( new CLHEP::RandGeneral( engine , &(binnedEnergySpectrum()[0]), _nbins) );
    } else {
      //  _shape = 0;
    }
  }

  void FoilParticleGenerator::generateFromFoil(ToyGenParticleCollection& genParts, long nParticles) {
    
    // Calculate spatial and temporal ranges
    _dcz  = (  _czmax -  _czmin);
    _dphi = ( _phimax - _phimin);
    _dt   = (   _tmax -   _tmin);
  

    // Get access to the geometry system.
    GeomHandle<Target> target;

    // Check if nfoils is bigger than 0;
    if (_nfoils < 1) {
      throw cms::Exception("GEOM")
        << "no foils are present";
    }

    for (int i=0; i<nParticles; i++) {
      
      // Pick a foil with a probability proportional to the volumes.
      int ifoil = static_cast<int>(_nfoils*_randFoils->fire());
      TargetFoil const& foil = target->foil(ifoil);
      
      // Foil properties.
      CLHEP::Hep3Vector const& center = foil.center();
      const double r1 = foil.rIn();
      const double dr = foil.rOut() - r1;
      
      // A random point within the foil. To change.
      const double r   = r1 + dr*_randFlat->fire();
      const double dz  = (-1.+2.*_randFlat->fire())*foil.halfThickness();
      const double phi = CLHEP::twopi*_randFlat->fire();
      CLHEP::Hep3Vector pos( center.x()+r*cos(phi), 
                             center.y()+r*sin(phi), 
                             center.z()+dz );
    
      // This should be an exponential decay.
      const double time = _tmin   +  _dt*(_randTime->fire());
      //_randFlat->fire();

      // Derived quantities.

      double e(0);
      CLHEP::Hep3Vector p3(0,0,0);
      CLHEP::HepLorentzVector mom(0,0,0,0);

      if (_genId == GenId::conversionGun) {
        e = sqrt( _p*_p + _mass*_mass );
      }

      if (_genId == GenId::dio1) {
        e  = _elow + _shape->fire() * (_ehi - _elow);
        _p = safeSqrt(e*e - _mass*_mass);
      }

      // Pick random 3 vector with the requested momentum.
      p3 = _randomUnitSphere->fire(_p);
      mom.setPx( p3.x() );
      mom.setPy( p3.y() );
      mom.setPz( p3.z() );
      mom.setE( e );
      
      // Add the particle to  the list.
      genParts.push_back( ToyGenParticle( _pdgId, _genId, pos, mom, time));

    } // end loop over generated DIO electrons
 
  }

  // Energy spectrum of the electron from DIO.
  double FoilParticleGenerator::energySpectrum( double e )
  {
    return pow<5>(_conversionEnergyAluminum - e) ;
  } 

  // Compute a binned representation of the energy spectrum of the electron from DIO.
  std::vector<double> FoilParticleGenerator::binnedEnergySpectrum(){

    // Sanity check.
    if (_nbins <= 0) {
      throw cms::Exception("RANGE") 
        << "Nonsense DecayInOrbitGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dE = (_ehi - _elow) / _nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);
    
    for (int ib=0; ib<_nbins; ib++) {
      double x = _elow+(ib+0.5) * dE;
      spectrum.push_back(energySpectrum(x));
    }

    return spectrum;
  } // FoilParticleGenerator::binnedEnergySpectrum

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
