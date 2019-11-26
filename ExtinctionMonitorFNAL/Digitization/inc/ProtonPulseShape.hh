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

#include "ExtinctionMonitorFNAL/Digitization/inc/PixelCharge.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"

namespace fhicl { class ParameterSet; }
namespace art   { class Event; }

namespace mu2e {

  class SimParticle;
  class SimParticleParentGetter;
  class MARSInfo;

  namespace ExtMonFNAL {

    class ProtonPulseShape {

      bool marsMode_;
      std::string marsInfoModuleLabel_;
      std::string marsInfoInstanceName_;

      // For pure G4 jobs hit times should be randomized based on the
      // "most primary" SimParticle (corresponding to a GenParticle)
      // in the chain .  However if we start with MARS inputs,
      // different GenParticles may correspond to the same primary
      // proton, so we need to go back to that proton.

      struct ProtonPathMARSId {
        MARSInfo minfo;
        IO::G4JobInfo g4s1info;
        ProtonPathMARSId(const MARSInfo& m, const IO::G4JobInfo& g) : minfo(m), g4s1info(g) {}
      };

      struct ProtonPathSort {
        bool operator()(const ProtonPathMARSId& a, const ProtonPathMARSId& b) const {
          CmpProtonIdAndSimPath cm;
          IO::CmpG4JobInfo cg;
          return cm(a.minfo,b.minfo) || (!cm(b.minfo, a.minfo) &&
                                         cg(a.g4s1info, b.g4s1info));
        }
      };

      typedef std::map<ProtonPathMARSId, double, ProtonPathSort> TimeMapMARS;
      TimeMapMARS tmm_;

      typedef std::map<art::Ptr<SimParticle>, double> TimeMapG4;
      TimeMapG4 tgm_;

      CLHEP::RandFlat flat_;
      double pulseHalfWidth_;
      bool messagePrinted_;

      void apply(PixelChargeHistory *inout, const SimParticleParentGetter& pg, const art::Event& event);

      double getTimeShiftForPrimary(const art::Ptr<SimParticle>& particle,
                                    const SimParticleParentGetter& pg,
                                    const art::Event& event);

      art::Ptr<SimParticle> getG4Primary(const art::Ptr<SimParticle>& particle,
                                         const SimParticleParentGetter& pg);

      ProtonPathMARSId getPrimaryMARSId(const art::Ptr<SimParticle>& particle, const art::Event& event);

    public:

      ProtonPulseShape(const fhicl::ParameterSet& pset,
                       art::RandomNumberGenerator::base_engine_t& rng);

      void apply(PixelChargeCollection *inout, const art::Event& event);
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Digitization_inc_ProtonPulseShape_hh*/
