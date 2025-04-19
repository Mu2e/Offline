//
// DefinesGeant4  scoring meshes and associated scorers for Mu2e geometry
//
// Original author BE
//

#include <cstdio>

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Principal/SubRun.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
//#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoringManager.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/MCDataProducts/inc/ScorerSummary.hh"

#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4PSPassageCellFlux.hh"
#include "G4ScoreLogColorMap.hh"
#include "G4PSCellFlux.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackCounter.hh"



namespace mu2e {

  Mu2eG4ScoringManager::Mu2eG4ScoringManager(G4ScoringManager* fSMan, const Mu2eG4Config::Scoring& config):
     fSMan_(fSMan),
     enabled_(config.enabled()),
     meshNames_(config.meshNames()),
     scorerNames_(config.scorerNames()),
     writeFile_(false),
     fileDirectory_("")
  {
    config.writeFile(writeFile_);
    config.fileDirectory(fileDirectory_);

    if (enabled_ && (meshNames_.empty() || scorerNames_.empty()))
       throw cet::exception("BADINPUT")<<"Mu2eG4ScoringManager: scoring requires at least one mesh name and one scorer name"
                                       << std::endl;
  }


  void Mu2eG4ScoringManager::initialize()
  {
    if (!enabled_) return;

    auto& geom   = *(art::ServiceHandle<GeometryService>());
    auto& config = geom.config();
    auto& reg    = (art::ServiceHandle<Mu2eG4Helper>())->antiLeakRegistry();

    const int verboseLevel = config.getInt("scoring.verboseLevel");
    const int nMesh        =  meshNames_.size();

    std::vector<double> meshPositionX, meshPositionY, meshPositionZ;
    config.getVectorDouble("scoring.meshPositionX", meshPositionX, nMesh);
    config.getVectorDouble("scoring.meshPositionY", meshPositionY, nMesh);
    config.getVectorDouble("scoring.meshPositionZ", meshPositionZ, nMesh);

    std::vector<double> meshSizeX, meshSizeY, meshSizeZ;
    config.getVectorDouble("scoring.meshSizeX", meshSizeX, nMesh);
    config.getVectorDouble("scoring.meshsizeY", meshSizeY, nMesh);
    config.getVectorDouble("scoring.meshsizeZ", meshSizeZ, nMesh);

    std::vector<int> meshSegmentX, meshSegmentY, meshSegmentZ;
    config.getVectorInt("scoring.meshSegmentX", meshSegmentX, nMesh);
    config.getVectorInt("scoring.meshSegmentY", meshSegmentY, nMesh);
    config.getVectorInt("scoring.meshSegmentZ", meshSegmentZ, nMesh);


    for (size_t i=0;i<meshNames_.size();++i){
      auto mesh = reg.add(new G4ScoringBox(meshNames_[i]));
      mesh->SetVerboseLevel(verboseLevel);

      G4double vsize[3] ={meshSizeX[i],meshSizeY[i],meshSizeZ[i]};
      mesh->SetSize(vsize);

      G4double centerPosition[3]={meshPositionX[i],meshPositionY[i],meshPositionZ[i]};
      mesh->SetCenterPosition(centerPosition);

      G4int nSegment[3]={meshSegmentX[i],meshSegmentY[i],meshSegmentZ[i]};
      mesh->SetNumberOfSegments(nSegment);

      for (const auto& psName: scorerNames_){
        switch (hashString(psName)) {
          case StringCode::CellFlux:
              mesh->SetPrimitiveScorer(new G4PSCellFlux(psName));
              break;
          case StringCode::FlatSurfaceFlux:
              mesh->SetPrimitiveScorer(new G4PSFlatSurfaceFlux(psName,0));
              break;
          case StringCode::DoseDeposit:
              mesh->SetPrimitiveScorer(new G4PSDoseDeposit(psName));
              break;
          case StringCode::EnergyDeposit:
              mesh->SetPrimitiveScorer(new G4PSEnergyDeposit(psName));
              break;
          case StringCode::TrackCounter:
              mesh->SetPrimitiveScorer(new G4PSTrackCounter(psName,1));
              break;
          default:
             throw cet::exception("BADINPUT")<<"Mu2eG4ScoringManager: unsupported scorer "<<psName
                                             <<". Choose among CellFlux, FlatSurfaceFlux, DoseDeposit, "
                                             <<"EnergyDeposit, TrackCounter\n"<< std::endl;
        }
      }

      fSMan_->RegisterScoringMesh(mesh);
      fSMan_->CloseCurrentMesh();
    }

    if (verboseLevel>1) fSMan_->List();
  }


  //------------------------------------------------------------------------------------------------------------
  void Mu2eG4ScoringManager::dumpInDataProduct(art::SubRun& subRun)
  {

    for (size_t i=0; i<fSMan_->GetNumberOfMesh(); ++i) {
      auto mesh     = fSMan_->GetMesh(i);
      auto meshName = mesh->GetWorldName();
      auto scoreMap = mesh->GetScoreMap();

      G4int nSegment[3]={0,0,0};
      mesh->GetNumberOfSegments(nSegment);

      G4VScoringMesh::MeshScoreMap::const_iterator msMapItr = scoreMap.begin();
      for (; msMapItr != scoreMap.end(); msMapItr++) {

        auto summaryColl = std::make_unique<ScorerSummaryCollection>();

        auto psName = msMapItr->first;
        auto score  = msMapItr->second->GetMap();
        auto unit   = mesh->GetPSUnitValue(psName);

        for (int ix = 0; ix < nSegment[0]; ++ix) {
          for (int iy = 0; iy < nSegment[1]; ++iy) {
            for (int iz = 0; iz < nSegment[2]; ++iz) {

              G4int idx      = ix*nSegment[1]*nSegment[2] +iy*nSegment[2]+iz; //check G4VScoreWriter for changes
              auto  value    = score->find(idx);
              int   entries  = (value != score->end()) ? value->second->n() : 0;
              float total    = (value != score->end()) ? value->second->sum_wx()/unit : 0.0;
              float totalSqr = (value != score->end()) ? value->second->sum_wx2()/unit : 0.0;

              summaryColl->emplace_back(ScorerSummary(ix,iy,iz,entries,total,totalSqr));
            }
          }
        }

        std::string instanceName = meshName + psName;
        subRun.put(std::move(summaryColl),instanceName,art::fullSubRun());
     }

     if (writeFile_) {
      std::string filename = fileDirectory_ + fSMan_->GetWorldName(i) + std::string(".dat");
      fSMan_->DumpAllQuantitiesToFile(fSMan_->GetWorldName(i),filename);
     }
   }

 }


  //------------------------------------------------------------------------------------------------------------
  void Mu2eG4ScoringManager::reset()
  {
    for (size_t i=0;i<fSMan_->GetNumberOfMesh();++i) fSMan_->GetMesh(i)->ResetScore();
  }


  //------------------------------------------------------------------------------------------------------------
  Mu2eG4ScoringManager::StringCode Mu2eG4ScoringManager::hashString(const G4String& str)
  {
    if (str == "CellFlux")        return StringCode::CellFlux;
    if (str == "FlatSurfaceFlux") return StringCode::FlatSurfaceFlux;
    if (str == "DoseDeposit")     return StringCode::DoseDeposit;
    if (str == "EnergyDeposit")   return StringCode::EnergyDeposit;
    if (str == "TrackCounter")    return StringCode::TrackCounter;
    return StringCode::Unknown;
  }

}
