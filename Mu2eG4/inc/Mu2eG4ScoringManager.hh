#ifndef Mu2eG4_ScoringManeger_hh
#define Mu2eG4_ScoringManager_hh
//
// Mu2eG4ScoringManager provides declarations for the Geant4 built-in scorer
// class for the Mu2e G4 simulation.
//
// Author: Bertrand Echenard
//

//G4 includes
#include "art/Framework/Principal/SubRun.h"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "G4ScoringManager.hh"
//#include "G4VScoreWriter.hh"


namespace mu2e {

  class Mu2eG4ScoringManager
  {
    public:
      enum class StringCode{CellFlux, FlatSurfaceFlux, DoseDeposit, EnergyDeposit,
                            TrackCounter, Unknown};

    public:
      Mu2eG4ScoringManager(G4ScoringManager* fSMan, const Mu2eG4Config::Scoring& config);
     ~Mu2eG4ScoringManager() = default;

      //methods
      void initialize();
      void dumpInDataProduct(art::SubRun& subRun);
      void reset();


    private:
      StringCode hashString(const G4String& str);

      G4ScoringManager*        fSMan_;     //non-owning G4 pointer
      bool                     enabled_;
      std::vector<std::string> meshNames_;
      std::vector<std::string> scorerNames_;
      bool                     writeFile_;
      std::string              fileDirectory_;
  };

}
#endif
