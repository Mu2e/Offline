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
      std::vector<elementData> elements;
      std::vector<double>      elementFractions;
    };
    
    GammaPairConversionSpectrum(){ _correlateAngleOverKE = false;}

    GammaPairConversionSpectrum(CLHEP::RandFlat* randFlat, bool correlateAngleOverKE = false);

    // random number generators are owned by the callers, no memory cleanup needed
    ~GammaPairConversionSpectrum(){}
   
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

  public:
    std::map<int, elementData> _elementMap; //map of element data by Z (standard A values assumed)
  private:

    CLHEP::RandFlat*   _rndFlat;
    bool               _correlateAngleOverKE; //add or remove correlation that disappeared in GEANT 4.10.4p03b --> 4.10.5.p01a
    int                _gMaxZet; //maximum element Z
    double             _me; // electron mass
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_GammaPairConversionSpectrum_hh */

