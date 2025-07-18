#ifndef Mu2eG4_ScoreWriter_hh
#define Mu2eG4_ScoreWriter_hh
//
// Mu2eG4ScoreWriter provides declarations to write the scores recorderded
// by the built-in Geant4 scorer
//
// Author: Bertrand Echenard
//


//G4 includes
#include "art/Framework/Principal/SubRun.h"
#include "G4VScoreWriter.hh"


namespace mu2e {

  class Mu2eG4ScoreWriter : public G4VScoreWriter
  {
    public:
      Mu2eG4ScoreWriter();
      virtual ~Mu2eG4ScoreWriter() = default;

      void dumpInDataProduct(art::SubRun& subRun);
      void dumpInFile(const std::string& fileDirectory);

    private:
      G4VScoringMesh* mesh_; //non-owning G4 pointer
  };

}
#endif
