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
#include "Offline/Mu2eG4/inc/Mu2eG4ScoreWriter.hh"
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

#include "CLHEP/Units/SystemOfUnits.h"


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

    std::vector<double> meshHalfSizeX, meshHalfSizeY, meshHalfSizeZ;
    config.getVectorDouble("scoring.meshHalfSizeX", meshHalfSizeX, nMesh);
    config.getVectorDouble("scoring.meshHalfSizeY", meshHalfSizeY, nMesh);
    config.getVectorDouble("scoring.meshHalfSizeZ", meshHalfSizeZ, nMesh);

    std::vector<int> meshSegmentX, meshSegmentY, meshSegmentZ;
    config.getVectorInt("scoring.meshSegmentX", meshSegmentX, nMesh);
    config.getVectorInt("scoring.meshSegmentY", meshSegmentY, nMesh);
    config.getVectorInt("scoring.meshSegmentZ", meshSegmentZ, nMesh);


    for (size_t i=0;i<meshNames_.size();++i){
      auto mesh = reg.add(new G4ScoringBox(meshNames_[i]));
      mesh->SetVerboseLevel(verboseLevel);

      G4double vsize[3] ={meshHalfSizeX[i]*CLHEP::mm,meshHalfSizeY[i]*CLHEP::mm,meshHalfSizeZ[i]*CLHEP::mm};
      mesh->SetSize(vsize);

      G4double centerPosition[3]={meshPositionX[i]*CLHEP::mm,meshPositionY[i]*CLHEP::mm,meshPositionZ[i]*CLHEP::mm};
      mesh->SetCenterPosition(centerPosition);

      G4int nSegment[3]={meshSegmentX[i],meshSegmentY[i],meshSegmentZ[i]};
      mesh->SetNumberOfSegments(nSegment);

      for (const auto& psName: scorerNames_){
        switch (hashString(psName)) {
          case StringCode::CellFlux:
              mesh->SetPrimitiveScorer(new G4PSCellFlux(psName));
              mesh->GetPrimitiveScorer(psName)->SetVerboseLevel(verboseLevel);
              break;
          case StringCode::FlatSurfaceFlux:
              mesh->SetPrimitiveScorer(new G4PSFlatSurfaceFlux(psName,0));
              mesh->GetPrimitiveScorer(psName)->SetVerboseLevel(verboseLevel);
              break;
          case StringCode::DoseDeposit:
              mesh->SetPrimitiveScorer(new G4PSDoseDeposit(psName));
              mesh->GetPrimitiveScorer(psName)->SetVerboseLevel(verboseLevel);
              break;
          case StringCode::EnergyDeposit:
              mesh->SetPrimitiveScorer(new G4PSEnergyDeposit(psName));
              mesh->GetPrimitiveScorer(psName)->SetVerboseLevel(verboseLevel);
              break;
          case StringCode::TrackCounter:
              mesh->SetPrimitiveScorer(new G4PSTrackCounter(psName,1));
              mesh->GetPrimitiveScorer(psName)->SetVerboseLevel(verboseLevel);
              break;
          default:
             throw cet::exception("BADINPUT")<<"Mu2eG4ScoringManager: unsupported scorer "<<psName
                                             <<". Choose among CellFlux, FlatSurfaceFlux, DoseDeposit, "
                                             <<"EnergyDeposit, TrackCounter\n"<< std::endl;
        }
      }
      mesh->SetVerboseLevel(verboseLevel);

      fSMan_->RegisterScoringMesh(mesh);
      fSMan_->CloseCurrentMesh();
    }

    if (verboseLevel>1) fSMan_->List();
  }


  //------------------------------------------------------------------------------------------------------------
  void Mu2eG4ScoringManager::dumpInDataProduct(art::SubRun& subRun)
  {
    for (size_t i=0; i<fSMan_->GetNumberOfMesh(); ++i) {
      Mu2eG4ScoreWriter writer(fSMan_->GetMesh(i));
      writer.dumpInDataProduct(subRun);
      if (writeFile_) writer.dumpInFile(fileDirectory_);
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
