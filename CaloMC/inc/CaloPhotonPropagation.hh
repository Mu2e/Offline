#ifndef CaloMC_CaloPhotonPropagation_hh
#define CaloMC_CaloPhotonPropagation_hh
//
// Sample the scintillation-photon propagation time from a location in the crystal.
// The time distribution vs. depth is taken from a detailed Geant4 simulation, stored
// as a 2D histogram (z on X, propagation time on Y) and turned into a per-depth CDF.
//
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include <string>
#include <vector>

namespace mu2e {

    class CaloPhotonPropagation
    {
        public:
          CaloPhotonPropagation(const std::string& fileName, const std::string& histName,
                                CLHEP::HepRandomEngine& engine);

          void  buildTable();
          float propTimeSimu(float z);        // sampled from the CDF (consumes a random number)
          float propTimeLine(float z) const;  // straight-line n*z/c estimate

        private:
          std::vector<float> timeProp_;      // Y-bin centers: candidate propagation times
          std::vector<float> cdf_;           // per-depth cumulative distribution, row-major [iz][itime]
          unsigned           nTimeDiv_{0};   // number of time bins  (histogram Y)
          unsigned           nZDiv_{0};      // number of depth bins (histogram X)
          float              dz_{0.f};       // depth-bin width
          CLHEP::RandFlat    randFlat_;
          std::string        fileName_;
          std::string        histName_;
          float              lightSpeed_{300.f};  // mm/ns, overwritten with c/n in buildTable
    };

}
#endif
