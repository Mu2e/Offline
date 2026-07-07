#ifndef CaloPhotonPropagation_HH
#define CaloPhotonPropagation_HH

// Calculate the propagation time from the location in the crystal
// Input based on detail Geant4 simulation of crystal

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include <vector>

namespace mu2e {

    class CaloPhotonPropagation
    {
        public:
          CaloPhotonPropagation(const std::string& fileName, const std::string& histName, CLHEP::HepRandomEngine& engine);

          void  buildTable  ();
          float propTimeSimu(float z);
          float propTimeLine(float z);

      private:
         std::vector<float>       timeProp_;
         std::vector<float>       cdf_;
         unsigned                 nTimeDiv_;
         unsigned                 nZDiv_;
         float                    dzTime_;
         CLHEP::RandFlat          randFlat_;
         std::string              fileName_;
         std::string              histName_;
         float                    lightSpeed_;
    };

}
#endif
