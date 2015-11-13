// Original author: Andrei Gaponenko, 2012

#ifndef Sources_inc_ExtMonFNALMARSUtils_hh
#define Sources_inc_ExtMonFNALMARSUtils_hh

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    // relevant information from a single line of the mars text file

    struct MARSParticle {
      unsigned protonNumber;
      int pid;
      double kineticEnergy;
      double weight;
      double x, y, z;
      double dcx, dcy, dcz;
      double tof;
      MARSParticle()
        : protonNumber(-1U), pid(), kineticEnergy(), weight(), x(), y(), z(), dcx(), dcy(), dcz(), tof()
      {}
    };

    std::ostream& operator<<(std::ostream& os, const MARSParticle& mp);

    bool readMARSLine(std::istream& file, MARSParticle& res);

    //================================================================
    class MARSMu2eConverter {
      GlobalConstantsHandle<ParticleDataTable> pdt_;
    public:
      GenParticle marsToMu2eParticle(const MARSParticle& mp);

      // returns pdgId
      PDGCode::type marsToMu2eParticleCode(int marsPID);

      inline CLHEP::Hep3Vector marsToMu2ePosition(double x, double y, double z) {
        // Apply the offsets and convert cm to mm
        return CLHEP::Hep3Vector(10*(x+390.4), 10*(y+0), 10*(z-906.8));
      }

      inline double marsToMu2eTime(double t) {
        // seconds to ns
        return t*1.e9;
      }

      inline double marsToMu2eEnergy(double ke) {
        // GeV to MeV
        return ke * 1.e3;
      }
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*Sources_inc_ExtMonFNALMARSUtils_hh*/
