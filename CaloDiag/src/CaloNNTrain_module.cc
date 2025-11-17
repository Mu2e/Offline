#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"

#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"



#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"

#include <vector>
#include <string>



namespace mu2e {

  class CaloNNTrain : public art::EDAnalyzer
  {
     public:
       struct Config
       {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Atom<art::InputTag> caloClusterCollection   { Name("caloClusterCollection"),   Comment("Calo cluster collection name") };
           fhicl::Atom<art::InputTag> caloClusterMCCollection { Name("caloClusterMCCollection"), Comment("Calo cluster MC collection name") };
           fhicl::Atom<float>         minEtoTest              { Name("minEtoTest"),              Comment("Minimum Energy to run the MVA") };
           fhicl::Atom<float>         MCEdepCut               { Name("MCEdepCut"),               Comment("Min E cut for MC contribution to cluster") };
           fhicl::Atom<int>           diagLevel               { Name("diagLevel"),               Comment("Diag Level"),0 };
       };

       explicit CaloNNTrain(const art::EDAnalyzer::Table<Config>& config) :
         EDAnalyzer{config},
         caloClusterToken_  {consumes<CaloClusterCollection> (config().caloClusterCollection())},
         caloClusterMCToken_{consumes<CaloClusterMCTruthAssn>(config().caloClusterMCCollection())},
         minEtoTest_        (config().minEtoTest()),
         MCEdepCut_         (config().MCEdepCut()),
         diagLevel_         (config().diagLevel())
       {}

       virtual void analyze(const art::Event& event) override;
       virtual void beginJob() override;


     private:
       float secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const;

       art::ProductToken<CaloClusterCollection>  caloClusterToken_;
       art::ProductToken<CaloClusterMCTruthAssn> caloClusterMCToken_;
       float   minEtoTest_;
       float   MCEdepCut_;
       int     diagLevel_;

       TTree* Ntup_;
       int    evt_,cluNcrys_,cluDisk_;
       float  cluEnergy_,cluTime_,cluCogR_,cluE1_,cluE2_,cluE9_,cluE25_,cluSec_,cluMCEtot_;
       std::vector<int> cluMCEPdg_,cluMCCr_;
       std::vector<float> cluMCEdep_;
  };


  void CaloNNTrain::beginJob()
  {
     art::ServiceHandle<art::TFileService> tfs;
     Ntup_  = tfs->make<TTree>("Calo", "Calo");
     Ntup_->Branch("evt",         &evt_ ,         "evt/I");
     Ntup_->Branch("cluNcrys",    &cluNcrys_ ,    "cluNcrys/I");
     Ntup_->Branch("cluDisk",     &cluDisk_ ,     "cluDisk/I");
     Ntup_->Branch("cluEnergy",   &cluEnergy_ ,   "cluEnergy/F");
     Ntup_->Branch("cluTime",     &cluTime_ ,     "cluTime/F");
     Ntup_->Branch("cluCogR",     &cluCogR_ ,     "cluCogR/F");
     Ntup_->Branch("cluE1",       &cluE1_ ,       "cluE1/F");
     Ntup_->Branch("cluE2",       &cluE2_ ,       "cluE2/F");
     Ntup_->Branch("cluE9",       &cluE9_ ,       "cluE9/F");
     Ntup_->Branch("cluE25",      &cluE25_ ,      "cluE25/F");
     Ntup_->Branch("cluSec",      &cluSec_ ,      "cluSec/F");
     Ntup_->Branch("cluMCEtot",   &cluMCEtot_ ,   "cluMCEtot/F");
     Ntup_->Branch("cluMCEdep",   &cluMCEdep_);
     Ntup_->Branch("cluMCEPdg",   &cluMCEPdg_);
     Ntup_->Branch("cluMCCr",     &cluMCCr_);
  }


  void CaloNNTrain::analyze(const art::Event& event)
  {
     evt_ = event.id().event();
     art::Handle<CaloClusterCollection>  caloClustersHandle   = event.getHandle<CaloClusterCollection>(caloClusterToken_);
     art::Handle<CaloClusterMCTruthAssn> caloClustersMCHandle = event.getHandle<CaloClusterMCTruthAssn>(caloClusterMCToken_);

     const Calorimeter& cal = *(GeomHandle<Calorimeter>());
     const CaloClusterCollection&  caloClusters(*caloClustersHandle);
     const CaloClusterMCTruthAssn& caloClustersMC(*caloClustersMCHandle);

     for (const auto& cluster : caloClusters)
     {
        if (cluster.energyDep() < minEtoTest_) continue;

        const auto& hits          = cluster.caloHitsPtrVector();
        const auto& neighborsId   = cal.crystal(hits[0]->crystalID()).neighbors();
        const auto& nneighborsId  = cal.crystal(hits[0]->crystalID()).nextNeighbors();


        float e1(hits[0]->energyDep()),e2(hits[0]->energyDep());
        if (hits.size()>1) e2 += hits[1]->energyDep();

        float e9(hits[0]->energyDep()),e25(hits[0]->energyDep());
        for (auto hit : hits)
        {
            if (std::find(neighborsId.begin(),  neighborsId.end(),  hit->crystalID()) != neighborsId.end())  {e9 += hit->energyDep();e25 += hit->energyDep();}
            if (std::find(nneighborsId.begin(), nneighborsId.end(), hit->crystalID()) != nneighborsId.end()) {e25 += hit->energyDep();}
        }

        float secondMom = secondMoment(cal,hits);

        //MC analyis
        float MCEdepTot(0);
        cluMCEdep_.clear();cluMCEPdg_.clear();cluMCCr_.clear();

        auto itMC = caloClustersMC.begin();
        while (itMC != caloClustersMC.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}
        const auto eDepMCs = (itMC != caloClustersMC.end()) ? itMC->second->energyDeposits() : std::vector<CaloEDepMC>{};

        if (itMC != caloClustersMC.end()){
          MCEdepTot= itMC->second->totalEnergyDep();
          for (auto& edep : itMC->second->energyDeposits()) {
             if (edep.energyDep() < MCEdepCut_) continue;
             cluMCEdep_.push_back(edep.energyDep());
             cluMCEPdg_.push_back(edep.sim()->pdgId());
             cluMCCr_.push_back(edep.sim()->creationCode());
          }
        }

       cluNcrys_  = cluster.size();
       cluDisk_   = cluster.diskID();
       cluEnergy_ = cluster.energyDep();
       cluTime_   = cluster.time();
       cluCogR_   = cluster.cog3Vector().perp();
       cluE1_     = e1/cluster.energyDep();
       cluE2_     = e2/cluster.energyDep();
       cluE9_     = e9/cluster.energyDep();
       cluE25_    = e25/cluster.energyDep();
       cluSec_    = secondMom;
       cluMCEtot_ = MCEdepTot;

       Ntup_->Fill();
    }
  }


  //-----------------------------------------------------------------------------------------------------------------------
  float CaloNNTrain::secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const
  {
     double sx(0),sy(0),sx2(0),sy2(0),sw(0);
     for (const auto& hit : hits){
        auto crId(hit->crystalID());
        auto energy(hit->energyDep());

        auto xCrystal = cal.crystal(crId).position().x();
        auto yCrystal = cal.crystal(crId).position().y();
        auto weight   = energy; //maybe try log(energy);

        sw  += weight;
        sx  += xCrystal*weight;
        sy  += yCrystal*weight;
        sx2 += xCrystal*xCrystal*weight;
        sy2 += yCrystal*yCrystal*weight;
     }
     return (sx2-sx*sx/sw + sy2-sy*sy/sw)/sw;
  }

}

DEFINE_ART_MODULE(mu2e::CaloNNTrain)
