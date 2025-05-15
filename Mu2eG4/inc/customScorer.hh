#ifndef customScorer_h
#define customScorer_h 1

//
// Custom scorer to debug Mu2eG4Scoring
// Author BE

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4PSDirectionFlag.hh"

namespace mu2e {

  class customScorer : public G4VPrimitiveScorer
  {
   public:
    customScorer(G4String name, G4int direction, G4int depth = 0);
    ~customScorer() override = default;

    inline void Weighted(G4bool flg = true) { weighted = flg; }

    void Initialize(G4HCofThisEvent*) override;
    void clear() override;
    void PrintAll() override;

    virtual void SetUnit(const G4String& unit);


   protected:
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    G4int GetIndex(G4Step*) override;


   private:
    void dumpStep(G4Step* aStep, G4VSolid* solid);

    G4int HCID{-1};
    G4int fDirection;
    G4THitsMap<G4double>* EvtMap{nullptr};
    G4bool weighted{false};
    G4int fDepthi, fDepthj, fDepthk;
  };

}
#endif
