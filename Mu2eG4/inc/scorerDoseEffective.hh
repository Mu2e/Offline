#ifndef scorerDoseEffective_h
#define scorerDoseEffective_h 1
//
// Custom scorer to calculate exposed dose using fluence-to-dose conversion tables
// Author BE
//
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "CLHEP/Random/RandFlat.h"

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4PSDirectionFlag.hh"


namespace mu2e {

  class scorerDoseEffective : public G4VPrimitiveScorer
  {
    public:
      scorerDoseEffective(const G4String& name,
                          const Mu2eG4Config::Physics& configPhysics,
                          G4int depth = 0);
      ~scorerDoseEffective() override = default;

      void Initialize(G4HCofThisEvent*) override;
      void clear()                      override;
      void PrintAll();

    protected:
      G4bool   ProcessHits(G4Step*, G4TouchableHistory*) override;
      G4int    GetIndex(G4Step*)                         override;
      G4double ComputeVolume(G4Step* aStep, G4int idx);

    private:
      G4int                 HCID_{-1};
      G4THitsMap<G4double>* EvtMap_{nullptr};
      G4int                 fDepthi_;
      G4int                 fDepthj_;
      G4int                 fDepthk_;
      scorerFTDConverter    FTDConverter_;
  };

}
#endif
