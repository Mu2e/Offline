// Mu2eG4ProcessCopier: see the header.

#include "Offline/Mu2eG4/inc/Mu2eG4ProcessCopier.hh"

#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4VParticleChange.hh"
#include "Geant4/G4VPhysicalVolume.hh"

namespace mu2e {

  Mu2eG4ProcessCopier::Mu2eG4ProcessCopier(G4VProcess* wrapped, unsigned elsewhere,
                                           std::vector<std::pair<const G4Region*, unsigned>> const& regions)
    : G4WrapperProcess("", wrapped->GetProcessType())   // RegisterProcess appends the wrapped name to ""
    , elsewhere_(elsewhere)
    , regions_(regions) {
    RegisterProcess(wrapped);
    SetProcessSubType(wrapped->GetProcessSubType());
    SetVerboseLevel(wrapped->GetVerboseLevel());
  }

  unsigned Mu2eG4ProcessCopier::factor(const G4Step& step) const {
    // a discrete interaction happens inside the volume of the pre-step point (a step that ends on a
    // boundary is limited by the geometry, not by this process)
    const G4VPhysicalVolume* pv = step.GetPreStepPoint()->GetPhysicalVolume();
    if (pv == nullptr) { return elsewhere_; }
    const G4Region* region = pv->GetLogicalVolume()->GetRegion();
    for (auto const& r : regions_) {
      if (r.first == region) { return r.second; }
    }
    return elsewhere_;
  }

  G4VParticleChange* Mu2eG4ProcessCopier::PostStepDoIt(const G4Track& track, const G4Step& step) {
    const unsigned n = factor(step);
    if (n <= 1) { return pRegProcess->PostStepDoIt(track, step); }

    // Each call re-initializes the wrapped process's particle change from the track, which Geant4 has not
    // updated yet, so every call samples the same interaction. Clear() only resets the count of
    // secondaries; the G4Track objects stay alive and are handed back below.
    saved_.clear();
    G4VParticleChange* pc = nullptr;
    for (unsigned i = 0; i < n; ++i) {
      pc = pRegProcess->PostStepDoIt(track, step);
      if (i + 1 < n) {
        for (G4int j = 0; j < pc->GetNumberOfSecondaries(); ++j) { saved_.push_back(pc->GetSecondary(j)); }
        pc->Clear();
      }
    }
    // pc holds the final state of the last call, with its own secondaries
    const double inv = 1. / static_cast<double>(n);
    for (G4int j = 0; j < pc->GetNumberOfSecondaries(); ++j) {
      G4Track* t = pc->GetSecondary(j);
      t->SetWeight(t->GetWeight() * inv);
    }
    const bool byProcess = pc->IsSecondaryWeightSetByProcess();
    pc->SetSecondaryWeightByProcess(true);   // otherwise AddSecondary sets the parent's weight
    for (G4Track* t : saved_) {
      t->SetWeight(t->GetWeight() * inv);
      pc->AddSecondary(t);
    }
    pc->SetSecondaryWeightByProcess(byProcess);
    saved_.clear();
    return pc;
  }

  void Mu2eG4ProcessCopier::PrepareWorkerPhysicsTable(const G4ParticleDefinition& particle) {
    pRegProcess->PrepareWorkerPhysicsTable(particle);
  }

  void Mu2eG4ProcessCopier::BuildWorkerPhysicsTable(const G4ParticleDefinition& particle) {
    pRegProcess->BuildWorkerPhysicsTable(particle);
  }

} // namespace mu2e
