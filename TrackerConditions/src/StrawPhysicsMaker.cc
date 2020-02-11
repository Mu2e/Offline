
#include "TrackerConditions/inc/StrawPhysicsMaker.hh"
#include <cmath>
#include <algorithm>


using namespace std;

namespace mu2e {

  StrawPhysics::ptr_t StrawPhysicsMaker::fromFcl(StrawDrift::cptr_t strawDrift) {
    
    // integrate the number of ionizations
    double ptot = 0.0;
    std::vector<double> nProb = _config.probPerCharge();
    std::vector<double> eionize = _config.ionizationEnergyTable();
    std::vector<double> EIonize(eionize.size());
    std::vector<double> intNProb(nProb.size());

    // renormalize the probs
    for(unsigned iprob=0;iprob< nProb.size(); ++iprob ) {
      ptot += nProb[iprob];
    }
    double norm = 1.0/ptot;
    double psum = 0.0;

    // integrate the ionization energy table
    const double MeV_per_eV = 1.0e-6;
    for(size_t iion=0; iion<eionize.size(); ++iion){
      double esum = 0.0;
      if(iion>0) esum = EIonize[iion-1];
      esum += MeV_per_eV*eionize[iion] + _config.ionizedElectronKE();// cumulative energy to free 'iion' electrons
      EIonize[iion] = esum;
    }

    // now compute the average charge and energy per ionization 
    // from these distributions
    double NAverage = 0.0;
    double EAverage = 0.0;

    for(unsigned iprob=0; iprob< nProb.size(); ++iprob ) {
      unsigned nele = iprob+1;
      double nprob = nProb[iprob]*norm;
      intNProb[iprob] = psum + nprob;
      psum += nprob;

      double ionizationEnergy = 0.0;
      if(nele >=1 && nele <= EIonize.size()) {
	ionizationEnergy = EIonize[nele-1];
      } else if(nele<1) {
	ionizationEnergy = 0.0;
      } else {
	ionizationEnergy = EIonize.back();
      }

      EAverage += nprob*ionizationEnergy/nele;
      NAverage += nprob*nele;
    }

    auto ptr = std::make_shared<StrawPhysics>(
		 EIonize,
		 _config.meanFreePath(), _config.ionizedElectronKE(), 
		 _config.electronCharge(),
		 intNProb,
		 _config.gasGain(), _config.polyaA(),
		 _config.gainRMSSlope(),_config.nGainGauss(),
		 _config.propagationVelocity(),
		 NAverage,
		 EAverage,
		 _config.clusterDriftPolynomial(),
		 _config.driftTimeVariance(),
		 _config.useNonLinearDrift(),
		 _config.bFieldOverride(),
		 strawDrift);
 

    return ptr;
  }
 
}


