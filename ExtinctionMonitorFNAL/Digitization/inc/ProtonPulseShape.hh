// G4 simulations may be run with proton time==const, then a pulse
// shape can be applied at the digitization stage using this class.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Digitization_inc_ProtonPulseShape_hh
#define ExtinctionMonitorFNAL_Digitization_inc_ProtonPulseShape_hh

#include <string>
#include <map>

#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/FindOne.h"

#include "ExtinctionMonitorFNAL/Digitization/inc/PixelCharge.hh"

namespace fhicl { class ParameterSet; }
namespace art   { class Event; }

namespace mu2e {

  class SimParticle;
  class SimParticleParentGetter;
  class MARSInfo;

  namespace ExtMonFNAL {

    class ProtonPulseShape {

      struct PrimaryMARSId {
        int run;
        int subrun;
        int proton;
        bool operator<(const PrimaryMARSId& b) const {
          return
            (run < b.run) || ((run==b.run) &&
                              ((subrun < b.subrun) || ((subrun == b.subrun)
                                                       &&(proton < b.proton))));
        }

        PrimaryMARSId(int r, int s, int p) : run(r), subrun(s), proton(p) {}
      };

      typedef std::map<PrimaryMARSId, double> TimeMapMARS;

      // For pure G4 jobs hit times should be randomized based on the
      // "most primary" SimParticle (corresponding to a GenParticle)
      // in the chain .  However if we start with MARS inputs,
      // different GenParticles may correspond to the same primary
      // proton, so we need to go back to that proton.

      bool marsMode_;
      std::string marsInfoModuleLabel_;
      std::string marsInfoInstanceName_;

      CLHEP::RandFlat flat_;
      double pulseHalfWidth_;
      bool messagePrinted_;

      TimeMapMARS tmm_;

      void apply(PixelChargeHistory *inout,
                 const SimParticleParentGetter& pg,
                 const art::FindOne<MARSInfo>& mif);

      double getTimeShiftForPrimary(const art::Ptr<SimParticle>& particle,
                                    const SimParticleParentGetter& pg,
                                    const art::FindOne<MARSInfo>& mif);

      PrimaryMARSId getPrimaryMARSId(const art::Ptr<SimParticle>& particle,
                                     const SimParticleParentGetter& pg,
                                     const art::FindOne<MARSInfo>& mif);

    public:

      ProtonPulseShape(const fhicl::ParameterSet& pset,
                       art::RandomNumberGenerator::base_engine_t& rng);

      void apply(PixelChargeCollection *inout, const art::Event& event);
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Digitization_inc_ProtonPulseShape_hh*/
