#ifndef scorerDelayedDose_h
#define scorerDelayedDose_h 1
//
// Custom scorer to calculate residual radiation dose using fluence-to-dose conversion tables
// Author BE
//
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4PSDirectionFlag.hh"


namespace mu2e {

  class scorerDelayedDose : public G4VPrimitiveScorer
  {
     public:
       scorerDelayedDose(const G4String& name,
                         const Mu2eG4Config::Physics& configPhysics,
                         G4int depth = 0);
       ~scorerDelayedDose() override = default;

       void Initialize(G4HCofThisEvent*) override;
       void clear()                      override;
       void PrintAll();


     protected:
       G4bool   ProcessHits(G4Step*, G4TouchableHistory*) override;
       G4int    GetIndex(G4Step*)                         override;
       G4double ComputeVolume(G4Step* aStep, G4int idx);


     private:
       void     readTimeProfile(std::string filename);
       G4bool   IsInDecayWindow(G4double time);

       G4int                 HCID_{-1};
       G4THitsMap<G4double>* EvtMap_{nullptr};
       G4int                 fDepthi_;
       G4int                 fDepthj_;
       G4int                 fDepthk_;
       scorerFTDConverter    FTDConverter_;
       G4bool                isBiased_{false};
       std::vector<G4double> bin_{};
       std::vector<G4double> profile_{};
       G4double              totalCoolTime_{0};
  };

}
#endif
