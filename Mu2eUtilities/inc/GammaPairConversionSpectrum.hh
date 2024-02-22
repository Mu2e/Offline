#ifndef Mu2eUtilities_GammaPairConversionSpectrum_hh
#define Mu2eUtilities_GammaPairConversionSpectrum_hh

// Mu2e includes
// #include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

// cetlib includes
#include "cetlib_except/exception.h"

// C++ includes
#include <cmath>
#include <vector>
#include <map>
#include <utility>


namespace CLHEP {
  class RandFlat;
}
namespace mu2e {

  class GammaPairConversionSpectrum {

  public:
    struct elementData {
      int    z;
      double z3; //z^(1/3)
      double logZ3; //log(z)/3
      double deltaMaxLow;
      double deltaMaxHigh;
      double fCoulomb;
    };
    struct materialData {
      std::vector<elementData> elementDatas;
      std::vector<double>      elementFractions;
    };

    GammaPairConversionSpectrum(CLHEP::RandFlat* randFlat, bool correlateAngleOverKE = false);

    void fire(const CLHEP::HepLorentzVector &photon, elementData &material,
       CLHEP::HepLorentzVector &electron, CLHEP::HepLorentzVector &positron);

    void fire(const CLHEP::HepLorentzVector &photon, materialData &material,  //will randomly select an element in the material
       CLHEP::HepLorentzVector &electron, CLHEP::HepLorentzVector &positron);

    void betheHeitlerModel(const CLHEP::HepLorentzVector &photon, elementData &material,
       CLHEP::HepLorentzVector &electron, CLHEP::HepLorentzVector &positron);

    double screenFunction1 (const double delta);
    double screenFunction2 (const double delta);
    void   screenFunction12(const double delta, double &f1, double  &f2);

    //angular distribution
    void samplePairDirections(const CLHEP::HepLorentzVector &photon, double electron_ke,
       double positron_ke, CLHEP::Hep3Vector &electron_dir,
       CLHEP::Hep3Vector &positron_dir);

    double sampleCosTheta(double ke);
    double sampleThetaU(double ke);

    void initializeElementData();
    std::map<int, elementData> GetElementMap(){return _elementMap;};

  private:
    std::map<int, elementData> _elementMap; //map of element data by Z (standard A values assumed)


    CLHEP::RandFlat*   _rndFlat;
    bool               _correlateAngleOverKE; //add or remove correlation that disappeared in GEANT 4.10.4p03b --> 4.10.5.p01a
    double             _me; // electron mass
    const double nele = 2.;
    const double eps0max = 0.5; // photon energy>2*me <=> me<= 0.5 * photonenergy
    const double min_gamma_energy = 2.; //MeV
    const double middle_energy = 50.; // MeV
    const double DF_const=136.; // delta factor constant
    const double FZ_const=8.; // FZ constant
    const double DM_const=4.; // delta min constant
    const double cubic_root=1./3.;
    const double deltamin=1.4;
    const double sf1a[3]={42.038,-8.29,0.958}; // screen 1 funtion pars for delta>deltamin
    const double sf1b[3]={42.184,7.444,-1.623}; // screen 1 funtion pars for delta<=deltamin
    const double sf2b[3]={41.326,5.848,-0.902}; // screen 1 funtion pars for delta<=deltamin
    const double u_const=2.;
    const int z_max=121;
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_GammaPairConversionSpectrum_hh */
