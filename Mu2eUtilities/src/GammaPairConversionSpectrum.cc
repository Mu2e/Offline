#include "Mu2eUtilities/inc/GammaPairConversionSpectrum.hh"
// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "cetlib/pow.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"


namespace mu2e {
  GammaPairConversionSpectrum::GammaPairConversionSpectrum(CLHEP::RandFlat* rndFlat, bool correlateAngleOverKE):
    _rndFlat             (rndFlat),
    _correlateAngleOverKE(correlateAngleOverKE),
    _gMaxZet             (120){

    GlobalConstantsHandle<ParticleDataTable> pdt;
    _me    = pdt->particle(PDGCode::e_minus ).ref().mass().value();
    //initialize some element data
    GammaPairConversionSpectrum::initializeElementData();
  }

  //get a random conversion event
  void GammaPairConversionSpectrum::fire(const CLHEP::HepLorentzVector &photon, 
					 GammaPairConversionSpectrum::elementData &material, 
					 CLHEP::HepLorentzVector &electron, 
					 CLHEP::HepLorentzVector &positron) {
    if(photon.e() < 2.*_me)
      throw cet::exception("ERROR") << "GammaPairConversion::" << __func__ 
				    << " Photon energy below conversion threshold!";
    if(material.z < 1)
      throw cet::exception("ERROR") << "GammaPairConversion::" << __func__ 
				    << " Conversion material has a Z < 1!";
    //sample the actual spectrum
    betheHeitlerModel(photon, material, electron, positron);
      
  }

  //select a random element and then get a random conversion event
  void GammaPairConversionSpectrum::fire(const CLHEP::HepLorentzVector &photon, 
					 GammaPairConversionSpectrum::materialData &material, 
					 CLHEP::HepLorentzVector &electron, 
					 CLHEP::HepLorentzVector &positron) {
    if(photon.e() < 2.*_me)
      throw cet::exception("ERROR") << "GammaPairConversion::" << __func__ 
				    << " Photon energy below conversion threshold!";
    if(material.elements.size() < 1)
      throw cet::exception("ERROR") << "GammaPairConversion::" << __func__ 
				    << " Conversion material has no elements!";
    //select a random element
    double r, frac;
    unsigned nElements = material.elements.size(), index;
    do {
      r = _rndFlat->fire();
      index = (unsigned) nElements*_rndFlat->fire();
      frac = material.elementFractions[index];
    } while(r > frac);
    //fire with the specific material
    fire(photon, material.elements[index], electron, positron);
      
  }

  //recommended spectrum to use for photon energies below 80 GeV
  //From G4BetheHeitlerModel::SampleSecondaries
  void GammaPairConversionSpectrum::betheHeitlerModel(const CLHEP::HepLorentzVector &photon,
						      GammaPairConversionSpectrum::elementData &material, 
						      CLHEP::HepLorentzVector &electron,
						      CLHEP::HepLorentzVector &positron){

    double photon_energy = photon.e();
    double eps0          = _me/photon_energy;
    //must have at least 2*_me energy to convert
    if(eps0 > 0.5) return;

    double eps;
    double small_energy = 2.; //MeV
    if(photon_energy < small_energy) {
      //simple case if low energy conversion
      eps = eps0 + (0.5-eps0)*_rndFlat->fire();
    } else {
      //higher energy conversion (but BetheHeitlerModel used for < 80 GeV)
      double middle_energy = 50.; // MeV
      double deltaFactor = 136.*eps0/(material.z3);
      double deltaMax = material.deltaMaxLow;
      double FZ = 8.*(material.logZ3);
      if(photon_energy > middle_energy) {
	FZ += 8.*(material.fCoulomb);
	deltaMax = material.deltaMaxHigh;
      }
      double deltaMin = 4.*deltaFactor;
    
      //limits of esp
      double epsp = 0.5 - 0.5*std::sqrt(1.-deltaMin/deltaMax);
      double epsMin = std::max(eps0, epsp);
      double epsRange = 0.5 - epsMin;

      //sample the energy rate (eps) of the created electron (or positron)
      double F10, F20;
      GammaPairConversionSpectrum::screenFunction12(deltaMin, F10, F20);
      F10 -= FZ;
      F20 -= FZ;
      const double normF1   = std::max(F10 * epsRange * epsRange, 0.);
      const double normF2   = std::max(1.5 * F20                , 0.);
      const double normCond = normF1 / (normF1 + normF2);
      //three random numbers per sampling
      double greject = 0.;
      double rndnum;
      do {
	rndnum = _rndFlat->fire();
	if(normCond > rndnum) {
	  rndnum = _rndFlat->fire();
	  eps = 0.5 - epsRange * pow(rndnum,1./3.); //G4Pow::A13(rand) = rand^(1/3)
	  const double delta = deltaFactor/(eps*(1.-eps));
	  greject = (GammaPairConversionSpectrum::screenFunction1(delta) - FZ)/F10;
	} else {
	  rndnum = _rndFlat->fire();
	  eps = epsMin + epsRange*rndnum;
	  const double delta = deltaFactor/(eps*(1.-eps));
	  greject = (GammaPairConversionSpectrum::screenFunction2(delta) - FZ)/F20;
	}
	rndnum = _rndFlat->fire();
      } while (greject < rndnum);

      //randomly select charges
      const int charge = (_rndFlat->fire() < 0.5) ? 1 : -1;

      //get kinematics
      const double positron_energy = (charge > 0) ? eps*photon_energy : (1.-eps)*photon_energy;
      const double electron_energy = photon_energy - positron_energy;
      const double positron_ke = std::max(0., positron_energy - _me);
      const double electron_ke = std::max(0., electron_energy - _me);
      //get direction of e+e-
      CLHEP::Hep3Vector positron_dir, electron_dir;
      //sample from anglular distribution from Modified Tsai spectrum
      GammaPairConversionSpectrum::samplePairDirections(photon, electron_ke, positron_ke,
							electron_dir, positron_dir);

      //set electron and positron lorentz vectors
      electron_dir.setMag(std::sqrt(electron_energy*electron_energy - _me*_me)); //use as p vector
      positron_dir.setMag(std::sqrt(positron_energy*positron_energy - _me*_me));
      positron.setE(positron_energy);
      positron.setVect(positron_dir);
      electron.setE(electron_energy);
      electron.setVect(electron_dir);
    }
  }

  //from G4BetheHeitlerModel::ScreenFunction1
  double GammaPairConversionSpectrum::screenFunction1(const double delta) {
    return (delta > 1.4) ? 42.038 - 8.29*log(delta + 0.958) : 42.184 - delta*(7.444 - 1.623*delta);
  }
  //from G4BetheHeitlerModel::ScreenFunction2
  double GammaPairConversionSpectrum::screenFunction2(const double delta) {
    return (delta > 1.4) ? 42.038 - 8.29*log(delta + 0.958) : 41.326 - delta*(5.848 - 0.902*delta);
  }
  //from G4BetheHeitlerModel::ScreenFunction12
  void GammaPairConversionSpectrum::screenFunction12(const double delta, double &f1, double &f2) {
    if (delta > 1.4) {
      f1 = 42.038 - 8.29*log(delta + 0.958);
      f2 = f1;
    } else {
      f1 = 42.184 - delta*(7.444 - 1.623*delta);
      f2 = 41.326 - delta*(5.848 - 0.902*delta); 
    }
  }

  //From G4ModifiedTsai::SamplePairDirections
  void GammaPairConversionSpectrum::samplePairDirections(const CLHEP::HepLorentzVector &photon, double electron_ke,
							 double positron_ke, CLHEP::Hep3Vector &electron_dir,
							 CLHEP::Hep3Vector &positron_dir) {

    CLHEP::Hep3Vector photon_dir(photon.vect());
    photon_dir.setMag(1.);

    double phi = CLHEP::twopi * _rndFlat->fire();
    double sinp = std::sin(phi);
    double cosp = std::cos(phi);

    double u(0.), theta(0.), cost(0.), sint(0.);
    if(_correlateAngleOverKE) {
      u = GammaPairConversionSpectrum::sampleThetaU(electron_ke);
      theta = u*_me/electron_ke;
      cost = cos(theta);
    } else
      cost = GammaPairConversionSpectrum::sampleCosTheta(electron_ke);
    sint = std::sqrt((1.-cost)*(1.+cost));

    electron_dir.set(sint*cosp, sint*sinp, cost);
    electron_dir.rotateUz(photon_dir);

    if(_correlateAngleOverKE) {
      theta = u*_me/positron_ke;
      cost = cos(theta);
    } else
      cost = GammaPairConversionSpectrum::sampleCosTheta(positron_ke);
    sint = std::sqrt((1.-cost)*(1.+cost));

    positron_dir.set(-sint*cosp, -sint*sinp, cost);
    positron_dir.rotateUz(photon_dir);
  }

  //From G4ModifiedTsai::SampleCosTheta(G4double kinEnergy)
  double GammaPairConversionSpectrum::sampleCosTheta(double ke) {
    double uMax = 2.*(1. + ke/_me);
    static const double a1 = 1.6;
    static const double a2 = a1/3.;
    static const double border = 0.25;
    double u(0.);
    do {
      double rndnum1 = _rndFlat->fire();
      const double rndnum2 = _rndFlat->fire();
      const double uu = -log(rndnum1*rndnum2);
      rndnum1 = _rndFlat->fire();
      u = (border > rndnum1) ? uu*a1 : uu*a2;
    } while(u > uMax);
    //in GEANT 4.10.p03b used cos(u*me/ke) ~ 1 - (u*me/ke)^2 / 2 = 1 - 2*u^2/ (2 ke/me)^2 ~ 1 - 2*u^2/(2*(1+ke/me))^2
    return 1. - 2.*u*u/uMax/uMax;
  }

  //From G4BetheHeitlerModel::SampleSecondaries, used in 4.10.4.p03b
  double GammaPairConversionSpectrum::sampleThetaU(double ke) {
    static const double a1 = 1.6;
    static const double a2 = a1/3.;
    static const double border = 0.25;
    // double uMax = M_PI*ke/_me; //if above this angle > pi
    double rndnum1 = _rndFlat->fire(), rndnum2 = _rndFlat->fire();
    const double uu = -log(rndnum1*rndnum2);
    rndnum1 = _rndFlat->fire();
    const double u = (border > rndnum1) ? uu*a1 : uu*a2;
    return u;
  }

  void GammaPairConversionSpectrum::initializeElementData() {
    for(int z = 1; z < 121; ++z) {
      double z3 = pow(z, 1./3.);
      //from G4Element::ComputeCoulombFactor()
      const double k1 = 0.0083, k2 = 0.20206, k3 = 0.0020, k4 = 0.0369;
      const double alpha = 1./137.036; //fine structure constant
      double az2 = (alpha*z)*(alpha*z); //should use effective z!
      double az4 = az2*az2;
      double fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
      double FZLow = 8.*log(z3);
      double FZHigh = FZLow + 8.*fCoulomb;
      GammaPairConversionSpectrum::elementData data;
      data.fCoulomb = fCoulomb;
      data.z = z;
      data.z3 = z3;
      data.logZ3 = log(z)/3.;
      data.deltaMaxLow  = exp((42.038 - FZLow )/8.29) - 0.958;
      data.deltaMaxHigh = exp((42.038 - FZHigh)/8.29) - 0.958;
      _elementMap[z] = data;
    }
  }
}
