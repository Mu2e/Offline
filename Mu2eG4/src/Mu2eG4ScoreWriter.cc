//
// Custom score writer to dump information in data product
//
// Original author BE
//

#include "art/Framework/Principal/SubRun.h"

#include "Offline/MCDataProducts/inc/ScorerSummary.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoreWriter.hh"

#include "G4VScoringMesh.hh"



namespace mu2e {

  Mu2eG4ScoreWriter::Mu2eG4ScoreWriter() {}


  //------------------------------------------------------------------------------------------------------------
  void Mu2eG4ScoreWriter::dumpInDataProduct(art::SubRun& subRun)
  {
      auto meshName = fScoringMesh->GetWorldName();
      auto scoreMap = fScoringMesh->GetScoreMap();

      G4VScoringMesh::MeshScoreMap::const_iterator msMapItr = scoreMap.begin();
      for (; msMapItr != scoreMap.end(); msMapItr++) {

        auto summaryColl = std::make_unique<ScorerSummaryCollection>();

        auto psName = msMapItr->first;
        auto score  = msMapItr->second->GetMap();
        auto unit   = fScoringMesh->GetPSUnitValue(psName);

        for (int ix = 0; ix < fNMeshSegments[0]; ix++) {
          for (int iy = 0; iy < fNMeshSegments[1]; iy++) {
            for (int iz = 0; iz < fNMeshSegments[2]; iz++) {

              G4int idx = GetIndex(ix, iy, iz);
              auto  value    = score->find(idx);
              int   entries  = (value != score->end()) ? value->second->n() : 0;
              float total    = (value != score->end()) ? value->second->sum_wx()/unit : 0.0;
              float totalSqr = (value != score->end()) ? value->second->sum_wx2()/unit/unit : 0.0;

              if (entries>0) summaryColl->emplace_back(ScorerSummary(ix,iy,iz,entries,total,totalSqr));
            }
          }
        }

        std::string instanceName = meshName + psName;
        subRun.put(std::move(summaryColl),instanceName,art::fullSubRun());
     }
   }

   void Mu2eG4ScoreWriter::dumpInFile(const std::string& fileDirectory){
     G4VScoreWriter writer;
     writer.SetScoringMesh(fScoringMesh);
     std::string fileName = fileDirectory + fScoringMesh->GetWorldName() + std::string(".dat");
     writer.DumpAllQuantitiesToFile(fileName, "");
   }

 }
