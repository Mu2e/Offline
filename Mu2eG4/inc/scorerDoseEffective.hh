#ifndef scorerDoseEffective_h
#define scorerDoseEffective_h 1
//
// Custom scorer to calculate exposed dose using fluence-to-dose conversion tables
// Author BE
//
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "CLHEP/Random/RandFlat.h"

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4PSDirectionFlag.hh"


namespace mu2e {

  class scorerDoseEffective : public G4VPrimitiveScorer
  {
     using vectorPairF = std::vector<std::pair<float,float>>;

     public:
       scorerDoseEffective(const G4String& name, G4int depth = 0);
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
