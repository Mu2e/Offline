//
// DefinesGeant4  scoring meshes and associated scorers for Mu2e geometry
//
// Original author BE
//

#include <cstdio>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Principal/SubRun.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
//#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/ScorerConfigSummary.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoringManager.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoreWriter.hh"
#include "Offline/Mu2eG4/inc/scorerDoseEffective.hh"
#include "Offline/Mu2eG4/inc/scorerDelayedDose.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"

#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4PSPassageCellFlux.hh"
#include "G4ScoreLogColorMap.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSTrackCounter3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSVolumeFlux3D.hh"
#include "G4ParticleTable.hh"




#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  Mu2eG4ScoringManager::Mu2eG4ScoringManager(G4ScoringManager* fSMan,
                                             const Mu2eG4Config::Scoring& configScore,
                                             const Mu2eG4Config::Physics& configPhysics):
     fSMan_(fSMan),
     configPhysics_(configPhysics),
     enabled_(configScore.enabled()),
     meshNames_(configScore.meshNames()),
     scorerNames_(configScore.scorerNames()),
     writeFile_(false),
     fileDirectory_("")
  {
    configScore.writeFile(writeFile_);
    configScore.fileDirectory(fileDirectory_);

    if (enabled_ && (meshNames_.empty() || scorerNames_.empty()))
       throw cet::exception("BADINPUT")<<"Mu2eG4ScoringManager: scoring requires at least one "<<
                                         "mesh name and one scorer name"<< std::endl;
  }


  void Mu2eG4ScoringManager::initialize()
  {
    if (!enabled_) return;

    auto& geom   = *(art::ServiceHandle<GeometryService>());
    auto& config = geom.config();
    auto& reg    = (art::ServiceHandle<Mu2eG4Helper>())->antiLeakRegistry();

    const int verboseLevel = config.getInt("scoring.verboseLevel");
    const auto nMesh        = meshNames_.size();

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

        //Select the scorer
        switch (hashScorer(psName)) {
          case ScorerCode::CellFlux:
              mesh->SetPrimitiveScorer(new G4PSCellFlux3D(psName));
              break;
          case ScorerCode::DoseDeposit:
              mesh->SetPrimitiveScorer(new G4PSDoseDeposit3D(psName));
              break;
          case ScorerCode::EnergyDeposit:
              mesh->SetPrimitiveScorer(new G4PSEnergyDeposit3D(psName));
              break;
          case ScorerCode::FlatSurfaceFlux:
              mesh->SetPrimitiveScorer(new G4PSFlatSurfaceFlux3D(psName,1));
              break;
          case ScorerCode::PassageCellFlux:
              mesh->SetPrimitiveScorer(new G4PSPassageCellFlux3D(psName,1));
              break;
          case ScorerCode::TrackCounter:
              mesh->SetPrimitiveScorer(new G4PSTrackCounter3D(psName,1));
              break;
          case ScorerCode::VolumeFlux:
              mesh->SetPrimitiveScorer(new G4PSVolumeFlux3D(psName,1));
              break;
          case ScorerCode::DoseEffective:
              mesh->SetPrimitiveScorer(new scorerDoseEffective(psName,configPhysics_,1));
              break;
          case ScorerCode::DelayedDose:
              mesh->SetPrimitiveScorer(new scorerDelayedDose(psName,configPhysics_,1));
              break;
          default:
             throw cet::exception("BADINPUT")<<"Mu2eG4ScoringManager: unsupported scorer "<<psName<<". "
                                             <<"Choose among CellFlux, DoseDeposit, EnergyDeposit, "
                                             <<"FlatSurfaceFlux, TrackCounter, PassageCellFlux, VolumeFlux, "
                                             <<"DoseEffective DelayedDose\n";
        }

        //optionaly add a particle filter
        switch (hashParticle(psName)) {
          case ParticleCode::Electron:
              mesh->SetFilter(new G4SDParticleFilter("Electron filter",std::vector<G4String>{"e-","e+"}));
              break;
          case ParticleCode::Photon:
              mesh->SetFilter(new G4SDParticleFilter("Photon filter", std::vector<G4String>{"gamma"}));
              break;
          case ParticleCode::Pion:
              mesh->SetFilter(new G4SDParticleFilter("Pion filter",std::vector<G4String>{"pi+","pi-"} ));
              break;
          case ParticleCode::Proton:
              mesh->SetFilter(new G4SDParticleFilter("Proton filter",std::vector<G4String>{"proton","anti_proton"} ));
              break;
          case ParticleCode::Neutron:
              mesh->SetFilter(new G4SDParticleFilter("Neutron filter",std::vector<G4String>{"neutron","anti_neutron"} ));
              break;
          default:
              break; //no filter applied by default
        }
      }

      fSMan_->RegisterScoringMesh(mesh);
      mesh->SetVerboseLevel(verboseLevel);
      fSMan_->CloseCurrentMesh();
    }

    if (verboseLevel>1) fSMan_->List();
  }


  //------------------------------------------------------------------------------------------------------------
  void Mu2eG4ScoringManager::dumpInDataProduct(art::SubRun& subRun)
  {

    if (!enabled_) return;

    //Write the configuration of each mesh
    auto summaryConfigColl = std::make_unique<ScorerConfigSummaryCollection>();
    for (size_t i=0; i<fSMan_->GetNumberOfMesh(); ++i) {
      G4int nBins[3];
      fSMan_->GetMesh(i)->GetNumberOfSegments(nBins);
      const auto& name   = fSMan_->GetMesh(i)->GetWorldName();
      const auto& size   = fSMan_->GetMesh(i)->GetSize();
      const auto& center = fSMan_->GetMesh(i)->GetTranslation();
      summaryConfigColl->emplace_back(ScorerConfigSummary(name, nBins[0],nBins[1],nBins[2],size,center));
    }
    subRun.put(std::move(summaryConfigColl),art::fullSubRun());

    //Write the content of each mesh
    for (size_t i=0; i<fSMan_->GetNumberOfMesh(); ++i) {
      Mu2eG4ScoreWriter writer;
      writer.SetScoringMesh(fSMan_->GetMesh(i));
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
  Mu2eG4ScoringManager::ScorerCode Mu2eG4ScoringManager::hashScorer(const G4String& str)
  {
    if (str.find("CellFlux")        != std::string::npos) return ScorerCode::CellFlux;
    if (str.find("FlatSurfaceFlux") != std::string::npos) return ScorerCode::FlatSurfaceFlux;
    if (str.find("DoseDeposit")     != std::string::npos) return ScorerCode::DoseDeposit;
    if (str.find("EnergyDeposit")   != std::string::npos) return ScorerCode::EnergyDeposit;
    if (str.find("TrackCounter")    != std::string::npos) return ScorerCode::TrackCounter;
    if (str.find("DoseEffective")   != std::string::npos) return ScorerCode::DoseEffective;
    if (str.find("DelayedDose")     != std::string::npos) return ScorerCode::DelayedDose;
    return ScorerCode::Unknown;
  }

  //------------------------------------------------------------------------------------------------------------
  Mu2eG4ScoringManager::ParticleCode Mu2eG4ScoringManager::hashParticle(const G4String& str)
  {
    if (str.find("Electron") != std::string::npos) return ParticleCode::Electron;
    if (str.find("Photon")   != std::string::npos) return ParticleCode::Photon;
    if (str.find("Pion")     != std::string::npos) return ParticleCode::Pion;
    if (str.find("Proton")   != std::string::npos) return ParticleCode::Proton;
    if (str.find("Neutron")  != std::string::npos) return ParticleCode::Neutron;
    return ParticleCode::Unknown;
  }

}
