#ifndef Mu2eG4_Mu2eG4ProcessCopier_hh
#define Mu2eG4_Mu2eG4ProcessCopier_hh
//
// Mu2eG4ProcessCopier wraps a Geant4 process and makes its discrete interaction (PostStepDoIt) N times
// where it happens. The secondaries of all N samples are kept, each with its weight divided by N; the
// primary continues with the final state of the last sample. Each sample is an independent draw of the
// same interaction, so every quantity carried by the secondaries keeps its expectation, and the primary
// is as after one unbiased interaction.
//
// Geant4 EM secondary biasing (G4EmBiasingManager, Mu2eG4 physics.secondaryBiasing) splits only
// interactions with exactly one secondary. This wrapper works for any process: pair production (two
// secondaries) or a decay in flight (the primary disappears). Decays at rest (AtRestDoIt) and the
// continuous part (AlongStepDoIt) are passed through unchanged.
//
// The wrapper keeps the name, type and subtype of the wrapped process, so the creation codes of the
// secondaries and the stopping code of the primary do not change. N depends on the G4Region of the volume
// in which the interaction happens; elsewhere the wrapped process runs once (N = 1).
//
// Configured by Mu2eG4 physics.processCopies (Mu2eG4CustomizationPhysicsConstructor).
//

#include <vector>

#include "Geant4/G4WrapperProcess.hh"

class G4Region;
class G4Track;

namespace mu2e {

  class Mu2eG4ProcessCopier : public G4WrapperProcess {
  public:
    // everywhere > 0: N in every volume (regions must be empty); otherwise N per G4Region, 1 elsewhere
    Mu2eG4ProcessCopier(G4VProcess* wrapped, unsigned everywhere,
                        std::vector<std::pair<const G4Region*, unsigned>> const& regions);

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step) override;

    // G4WrapperProcess forwards these to the plain (master) versions of the wrapped process
    void PrepareWorkerPhysicsTable(const G4ParticleDefinition& particle) override;
    void BuildWorkerPhysicsTable(const G4ParticleDefinition& particle) override;

    unsigned factor(const G4Step& step) const;

  private:
    unsigned everywhere_;
    std::vector<std::pair<const G4Region*, unsigned>> regions_;
    std::vector<G4Track*> saved_;
  };

} // namespace mu2e

#endif
