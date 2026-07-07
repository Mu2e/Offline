//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"


namespace
{
   constexpr int ntupLen = 16384;
}


namespace mu2e {

  class CaloNNDiag : public art::EDAnalyzer {

     public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Atom<art::InputTag>     vdCollection          { Name("vdCollection"),           Comment("Virtual detector collection name") };
             fhicl::Atom<art::InputTag>     caloHitCollection     { Name("caloHitCollection"),      Comment("Calo Hit collection name") };
             fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"),  Comment("Calo cluster collection name") };
             fhicl::Atom<art::InputTag>     caloHitTruth          { Name("caloHitTruth"),           Comment("CaloHit truth name") };
             fhicl::Atom<art::InputTag>     caloClusterTruth      { Name("caloClusterTruth"),       Comment("caloCluster truth name") };
             fhicl::Atom<float>             minCluEnergy          { Name("minCluEnergy"),           Comment("Min cluster energy to write cluster") };
             fhicl::Atom<int>               diagLevel             { Name("diagLevel"),              Comment("Diag Level"),0 };
         };

       explicit CaloNNDiag(const art::EDAnalyzer::Table<Config>& config);
       virtual ~CaloNNDiag() {}

       virtual void beginJob();
       virtual void endJob() {};
       virtual void analyze(const art::Event& e);


     private:
       art::InputTag  virtualDetectorTag_;
       art::InputTag  caloHitTag_;
       art::InputTag  caloClusterTag_;
       art::InputTag  caloHitTruthTag_;
       art::InputTag  caloClusterTruthTag_;
       float          minCluEnergy_;
       int            diagLevel_;
       int            nProcess_;


       TTree* Ntup_;
       int   _evt,_run;

       int   nHits_,cryId_[ntupLen],crySectionId_[ntupLen],crySimIdx_[ntupLen],crySimLen_[ntupLen],cryConv_[ntupLen];
       float cryTime_[ntupLen],cryEdep_[ntupLen],cryEdepErr_[ntupLen],cryPosX_[ntupLen],cryPosY_[ntupLen],cryPosZ_[ntupLen];
       float cryTimeErr_[ntupLen],cryT1_[ntupLen],cryT2_[ntupLen],cryT1Err_[ntupLen],cryT2Err_[ntupLen];

       std::vector<std::vector<int>>   crySimId_,crySimPdgId_,crySimCrCode_;
       std::vector<std::vector<float>> crySimMom_,crySimTime_,crySimEdep_,crySimStartX_,crySimStartY_,crySimStartZ_;

       int   nCluster_,nCluSim_,cluNcrys_[ntupLen],cluDisk_[ntupLen];
       float cluEnergy_[ntupLen],cluEnergyErr_[ntupLen],cluTime_[ntupLen],cluTimeErr_[ntupLen],cluCogX_[ntupLen],cluCogY_[ntupLen],cluCogR_[ntupLen],
             cluCogZ_[ntupLen],cluE1_[ntupLen],cluE2_[ntupLen],cluE9_[ntupLen],cluE25_[ntupLen],cluSecMom_[ntupLen],cluDR_[ntupLen];
       int   cluSplit_[ntupLen],cluConv_[ntupLen],cluSimIdx_[ntupLen],cluSimLen_[ntupLen];

       std::vector<std::vector<int> > cluList_;
       std::vector<std::vector<int>>   cluSimId_,cluSimPdgId_,cluSimCrCode_;
       std::vector<std::vector<float>> cluSimMom_,cluSimMom2_,cluSimTime_,cluSimEdep_,cluSimPosX_,cluSimPosY_,cluSimPosZ_;
       std::vector<std::vector<float>> cluSimStartX_,cluSimStartY_,cluSimStartZ_;

       int   nVd_,vdId_[ntupLen],vdPdgId_[ntupLen],vdGenId_[ntupLen],vdGenIdx_[ntupLen];
       float vdTime_[ntupLen],vdPosX_[ntupLen],vdPosY_[ntupLen],vdPosZ_[ntupLen],vdMom_[ntupLen],vdMomX_[ntupLen],vdMomY_[ntupLen],vdMomZ_[ntupLen];
  };


  CaloNNDiag::CaloNNDiag(const art::EDAnalyzer::Table<Config>& config) :
    EDAnalyzer{config},
    virtualDetectorTag_ (config().vdCollection()),
    caloHitTag_         (config().caloHitCollection()),
    caloClusterTag_     (config().caloClusterCollection()),
    caloHitTruthTag_    (config().caloHitTruth()),
    caloClusterTruthTag_(config().caloClusterTruth()),
    minCluEnergy_       (config().minCluEnergy()),
    diagLevel_          (config().diagLevel()),
    nProcess_(0),
    Ntup_(0)
  {}

  void CaloNNDiag::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       Ntup_  = tfs->make<TTree>("Calo", "Calo");

       Ntup_->Branch("evt",          &_evt ,         "evt/I");
       Ntup_->Branch("run",          &_run ,         "run/I");

       Ntup_->Branch("nCry",         &nHits_ ,       "nCry/I");
       Ntup_->Branch("cryId",        &cryId_ ,       "cryId[nCry]/I");
       Ntup_->Branch("crySectionId", &crySectionId_, "crySectionId[nCry]/I");
       Ntup_->Branch("cryPosX",      &cryPosX_ ,     "cryPosX[nCry]/F");
       Ntup_->Branch("cryPosY",      &cryPosY_ ,     "cryPosY[nCry]/F");
       Ntup_->Branch("cryPosZ",      &cryPosZ_ ,     "cryPosZ[nCry]/F");
       Ntup_->Branch("cryEdep",      &cryEdep_ ,     "cryEdep[nCry]/F");
       Ntup_->Branch("cryEdepErr",   &cryEdepErr_ ,  "cryEdepErr[nCry]/F");
       Ntup_->Branch("cryTime",      &cryTime_ ,     "cryTime[nCry]/F");
       Ntup_->Branch("cryTimeErr",   &cryTimeErr_ ,  "cryTimeErr[nCry]/F");
       Ntup_->Branch("cryT1",        &cryT1_ ,       "cryT1[nCry]/F");
       Ntup_->Branch("cryT2",        &cryT2_ ,       "cryT2[nCry]/F");
       Ntup_->Branch("cryT1Err",     &cryT1Err_ ,    "cryT1Err[nCry]/F");
       Ntup_->Branch("cryT2Err",     &cryT2Err_ ,    "cryT2Err[nCry]/F");
       Ntup_->Branch("cryConv",      &cryConv_ ,     "cryConv[nCry]/I");
       Ntup_->Branch("crySimId",     &crySimId_);
       Ntup_->Branch("crySimPdgId",  &crySimPdgId_);
       Ntup_->Branch("crySimCrCode", &crySimCrCode_);
       Ntup_->Branch("crySimMom",    &crySimMom_);
       Ntup_->Branch("crySimTime",   &crySimTime_);
       Ntup_->Branch("crySimEdep",   &crySimEdep_);
       Ntup_->Branch("crySimStartX", &crySimStartX_);
       Ntup_->Branch("crySimStartY", &crySimStartY_);
       Ntup_->Branch("crySimStartZ", &crySimStartZ_);

       Ntup_->Branch("nCluster",     &nCluster_ ,    "nCluster/I");
       Ntup_->Branch("cluEnergy",    &cluEnergy_ ,   "cluEnergy[nCluster]/F");
       Ntup_->Branch("cluDisk",      &cluDisk_ ,     "cluDisk[nCluster]/I");
       Ntup_->Branch("cluEnergyErr", &cluEnergyErr_ ,"cluEnergyErr[nCluster]/F");
       Ntup_->Branch("cluTime",      &cluTime_ ,     "cluTime[nCluster]/F");
       Ntup_->Branch("cluTimeErr",   &cluTimeErr_ ,  "cluTimeErr[nCluster]/F");
       Ntup_->Branch("cluCogX",      &cluCogX_ ,     "cluCogX[nCluster]/F");
       Ntup_->Branch("cluCogY",      &cluCogY_ ,     "cluCogY[nCluster]/F");
       Ntup_->Branch("cluCogR",      &cluCogR_ ,     "cluCogR[nCluster]/F");
       Ntup_->Branch("cluCogZ",      &cluCogZ_ ,     "cluCogZ[nCluster]/F");
       Ntup_->Branch("cluNcrys",     &cluNcrys_ ,    "cluNcrys[nCluster]/I");
       Ntup_->Branch("cluE1",        &cluE1_ ,       "cluE1[nCluster]/F");
       Ntup_->Branch("cluE2",        &cluE2_ ,       "cluE2[nCluster]/F");
       Ntup_->Branch("cluE9",        &cluE9_ ,       "cluE9[nCluster]/F");
       Ntup_->Branch("cluE25",       &cluE25_ ,      "cluE25[nCluster]/F");
       Ntup_->Branch("cluSecMom",    &cluSecMom_ ,   "cluSecMom[nCluster]/F");
       Ntup_->Branch("cluSplit",     &cluSplit_ ,    "cluSplit[nCluster]/I");
       Ntup_->Branch("cluConv",      &cluConv_ ,     "cluConv[nCluster]/I");
       Ntup_->Branch("cluList",      &cluList_);

       Ntup_->Branch("cluSimId",     &cluSimId_);
       Ntup_->Branch("cluSimPdgId",  &cluSimPdgId_);
       Ntup_->Branch("cluSimCrCode", &cluSimCrCode_);
       Ntup_->Branch("cluSimMom",    &cluSimMom_);
       Ntup_->Branch("cluSimMom2",   &cluSimMom2_);
       Ntup_->Branch("cluSimTime",   &cluSimTime_);
       Ntup_->Branch("cluSimEdep",   &cluSimEdep_);
       Ntup_->Branch("cluSimPosX",   &cluSimPosX_);
       Ntup_->Branch("cluSimPosY",   &cluSimPosY_);
       Ntup_->Branch("cluSimPosZ",   &cluSimPosZ_);
       Ntup_->Branch("cluSimStartX", &cluSimStartX_);
       Ntup_->Branch("cluSimStartY", &cluSimStartY_);
       Ntup_->Branch("cluSimStartZ", &cluSimStartZ_);

       Ntup_->Branch("nVd",          &nVd_ ,         "nVd/I");
       Ntup_->Branch("vdId",         &vdId_ ,        "vdId[nVd]/I");
       Ntup_->Branch("vdPdgId",      &vdPdgId_ ,     "vdPdgId[nVd]/I");
       Ntup_->Branch("vdGenId",      &vdGenId_ ,     "vdGenId[nVd]/I");
       Ntup_->Branch("vdMom",        &vdMom_ ,       "vdMom[nVd]/F");
       Ntup_->Branch("vdMomX",       &vdMomX_ ,      "vdMomX[nVd]/F");
       Ntup_->Branch("vdMomY",       &vdMomY_ ,      "vdMomY[nVd]/F");
       Ntup_->Branch("vdMomZ",       &vdMomZ_ ,      "vdMomZ[nVd]/F");
       Ntup_->Branch("vdPosX",       &vdPosX_ ,      "vdPosX[nVd]/F");
       Ntup_->Branch("vdPosY",       &vdPosY_ ,      "vdPosY[nVd]/F");
       Ntup_->Branch("vdPosZ",       &vdPosZ_ ,      "vdPosZ[nVd]/F");
       Ntup_->Branch("vdTime",       &vdTime_ ,      "vdTime[nVd]/F");
       Ntup_->Branch("vdGenIdx",     &vdGenIdx_ ,    "vdGenIdx[nVd]/I");
  }



  void CaloNNDiag::analyze(const art::Event& event)
  {
      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloNNDiag =  "<<nProcess_ <<std::endl;

      double _mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if (!geom->hasElement<Calorimeter>() ) return;
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(caloHitTag_, CaloHitsHandle);
      if (!CaloHitsHandle.isValid()) return;
      const CaloHitCollection& CaloHits(*CaloHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(caloClusterTag_, caloClustersHandle);
      if (!caloClustersHandle.isValid()) return;
      const CaloClusterCollection& caloClusters(*caloClustersHandle);

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(virtualDetectorTag_,vdhits);

      //Calo digi truth assignment
      art::Handle<CaloHitMCTruthAssn> caloDigiTruthHandle;
      event.getByLabel(caloHitTruthTag_, caloDigiTruthHandle);
      if (!caloDigiTruthHandle.isValid()) return;
      const CaloHitMCTruthAssn& caloDigiTruth(*caloDigiTruthHandle);

      //Calo cluster truth assignment
      art::Handle<CaloClusterMCTruthAssn> caloClusterTruthHandle;
      event.getByLabel(caloClusterTruthTag_, caloClusterTruthHandle);
      if (!caloClusterTruthHandle.isValid()) return;
      const CaloClusterMCTruthAssn& caloClusterTruth(*caloClusterTruthHandle);


      //SELECT EVENTS WITH AT LEAST A CLUSTER ABOVE MINCLUENERGY_
      bool select(false);
      for (const auto& cluster : caloClusters) {
         if (cluster.energyDep()>minCluEnergy_) {select=true; break;}
      }
      if (!select) return;



      std::map<art::Ptr<SimParticle>, const StepPointMC*> vdMap;
      if (vdhits.isValid())
      {
         for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
         {
            const StepPointMC& hit = *iter;
            if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;
            vdMap[hit.simParticle()] = &hit;
         }
      }


       //--------------------------  Start --------------------------------
       _evt = event.id().event();
       _run = event.run();

       if (diagLevel_ == 3){std::cout << "processing event in calo_example " << nProcess_ << " run and event  = " << _run << " " << _evt << std::endl;}


       //--------------------------  Do calorimeter hits --------------------------------

       nHits_ = 0;
       crySimId_.clear();crySimPdgId_.clear();crySimCrCode_.clear();crySimTime_.clear();crySimEdep_.clear();
       crySimMom_.clear();crySimStartX_.clear();crySimStartY_.clear();crySimStartZ_.clear();

       for (unsigned int ic=0; ic<CaloHits.size();++ic){
           const CaloHit& hit            = CaloHits.at(ic);
           int diskId                    = cal.crystal(hit.crystalID()).diskID();
           CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.crystalID()).position());  //in disk FF frame

           //Find the caloDigiMC in the truth map
           auto itMC = caloDigiTruth.begin();
           while (itMC != caloDigiTruth.end()) {if (itMC->first.get() == &hit) break; ++itMC;}
           unsigned nCrySims = (itMC != caloDigiTruth.end()) ? itMC->second->nParticles() : 0;

           bool isConversion(false);
           if (nCrySims>0)
           {
              for (auto& edep : itMC->second->energyDeposits())
              {
                 auto parent(edep.sim());
                 while (parent->hasParent()) parent = parent->parent();
                 if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) isConversion=true;
              }
           }

           constexpr float invalid(999.0);
           float cryT1(invalid),cryT2(invalid),cryT1Err(invalid),cryT2Err(invalid);
           if (hit.recoCaloDigis().size()>1)
           {
              int idx0 = CaloSiPMId(hit.recoCaloDigis().at(0)->SiPMID()).SiPMLocalId();
              int idx1 = CaloSiPMId(hit.recoCaloDigis().at(1)->SiPMID()).SiPMLocalId();
              cryT1    = hit.recoCaloDigis().at(idx0)->time();
              cryT2    = hit.recoCaloDigis().at(idx1)->time();
              cryT1Err = hit.recoCaloDigis().at(idx0)->timeErr();
              cryT2Err = hit.recoCaloDigis().at(idx1)->timeErr();
           }

           cryId_[nHits_]        = hit.crystalID();
           crySectionId_[nHits_] = diskId;
           cryEdep_[nHits_]      = hit.energyDep();
           cryEdepErr_[nHits_]   = hit.energyDepErr();
           cryTime_[nHits_]      = hit.time();
           cryTimeErr_[nHits_]   = hit.timeErr();
           cryT1_[nHits_]        = cryT1;
           cryT2_[nHits_]        = cryT2;
           cryT1Err_[nHits_]     = cryT1Err;
           cryT2Err_[nHits_]     = cryT2Err;
           cryPosX_[nHits_]      = crystalPos.x();
           cryPosY_[nHits_]      = crystalPos.y();
           cryPosZ_[nHits_]      = crystalPos.z();
           cryConv_[nHits_]      = isConversion ? 1 : 0;

           std::vector<int> sId,sPdg,sCr;
           std::vector<float>st,sed,sm,stx,sty,stz;
           for (unsigned i=0;i< nCrySims;++i)
           {
               const auto& eDepMC = itMC->second->energyDeposit(i);

               auto parent(eDepMC.sim());
               while (parent->hasParent()) parent = parent->parent();

               sId.push_back(eDepMC.sim()->id().asInt());
               sPdg.push_back(eDepMC.sim()->pdgId());
               sCr.push_back(eDepMC.sim()->creationCode());
               sm.push_back(eDepMC.momentumIn());
               sed.push_back(eDepMC.energyDep());
               st.push_back(eDepMC.time());
               stx.push_back(parent->startPosition().x());
               sty.push_back(parent->startPosition().y());
               stz.push_back(parent->startPosition().z());
           }
           crySimId_.push_back(sId);
           crySimPdgId_.push_back(sPdg);
           crySimCrCode_.push_back(sCr);
           crySimMom_.push_back(sm);
           crySimEdep_.push_back(sed);
           crySimTime_.push_back(st);
           crySimStartX_.push_back(stx);
           crySimStartY_.push_back(sty);
           crySimStartZ_.push_back(stz);

           ++nHits_;
       }

       //--------------------------  Do clusters --------------------------------
       nCluster_ = 0;
       cluList_.clear();
       cluSimId_.clear();cluSimPdgId_.clear();cluSimCrCode_.clear();cluSimTime_.clear();cluSimEdep_.clear();
       cluSimMom_.clear();cluSimMom2_.clear();cluSimPosX_.clear();cluSimPosY_.clear();cluSimPosZ_.clear();
       cluSimStartX_.clear();cluSimStartY_.clear();cluSimStartZ_.clear();

       for (unsigned  ic=0; ic<caloClusters.size();++ic)
       {
          const CaloCluster& cluster = caloClusters.at(ic);
          if (cluster.energyDep() < minCluEnergy_) continue;

          std::vector<int> cryList;
          for (auto cryPtr : cluster.caloHitsPtrVector()) cryList.push_back(std::distance(&CaloHits.at(0),cryPtr.get()));

          ClusterUtils cluUtil(cal, cluster);

          auto itMC = caloClusterTruth.begin();
          while (itMC != caloClusterTruth.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}
          const auto eDepMCs = (itMC != caloClusterTruth.end()) ? itMC->second->energyDeposits() : std::vector<CaloEDepMC>{};

          bool isConv(false);
          if (itMC != caloClusterTruth.end()){
             for (auto& edep : itMC->second->energyDeposits()){
                if (edep.sim()->creationCode() == ProcessCode::mu2eCeMinusEndpoint ||
                    edep.sim()->creationCode() == ProcessCode::mu2eCePlusEndpoint ||
                    edep.sim()->creationCode() == ProcessCode::mu2eCeMinusLeadingLog ||
                    edep.sim()->creationCode() == ProcessCode::mu2eCePlusLeadingLog){
                  isConv = true;
                  break;
                }
             }
          }

          const CaloHit& seedHit    = CaloHits.at(cryList[0]);
          CLHEP::Hep3Vector seedPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(),cal.crystal(seedHit.crystalID()).position());

          cluDisk_[nCluster_]      = cluster.diskID();
          cluEnergy_[nCluster_]    = cluster.energyDep();
          cluEnergyErr_[nCluster_] = cluster.energyDepErr();
          cluTime_[nCluster_]      = cluster.time();
          cluTimeErr_[nCluster_]   = cluster.timeErr();
          cluNcrys_[nCluster_]     = cluster.size();
          cluCogX_[nCluster_]      = cluster.cog3Vector().x(); //in disk FF frame
          cluCogY_[nCluster_]      = cluster.cog3Vector().y();
          cluCogR_[nCluster_]      = sqrt(cluCogX_[nCluster_]*cluCogX_[nCluster_]+cluCogY_[nCluster_]*cluCogY_[nCluster_]);
          cluCogZ_[nCluster_]      = cluster.cog3Vector().z();
          cluE1_[nCluster_]        = cluUtil.e1();
          cluE2_[nCluster_]        = cluUtil.e2();
          cluE9_[nCluster_]        = cluUtil.e9();
          cluE25_[nCluster_]       = cluUtil.e25();
          cluSecMom_[nCluster_]    = cluUtil.secondMoment();
          cluSplit_[nCluster_]     = cluster.isSplit();
          cluConv_[nCluster_]      = isConv;//ic == icMCIdx;
          cluList_.push_back(cryList);


          std::vector<int> sId,sPdg,sCr;
          std::vector<float>st,sed,sm,sm2,spx,spy,spz,stx,sty,stz;
          for (unsigned i=0;i< eDepMCs.size();++i)
          {
              const auto& eDepMC = eDepMCs[i];
              art::Ptr<SimParticle> sim = eDepMC.sim();

              art::Ptr<SimParticle> smother(sim);
              while (smother->hasParent() && !smother->genParticle() ) smother = smother->parent();

              double simMom(-1);
              CLHEP::Hep3Vector simPos(0,0,0);
              auto vdMapEntry = vdMap.find(sim);
              if (vdMapEntry != vdMap.end())
              {
                 simMom = vdMapEntry->second->momentum().mag();
                 simPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(), vdMapEntry->second->position());
              }

              sId.push_back(sim->id().asInt());
              sPdg.push_back(sim->pdgId());
              sCr.push_back(sim->creationCode());
              st.push_back(eDepMC.time());
              sed.push_back(eDepMC.energyDep());
              sm.push_back(eDepMC.momentumIn());
              sm2.push_back(simMom);
              spx.push_back(simPos.x()); // in disk FF frame
              spy.push_back(simPos.y());
              spz.push_back(simPos.z());
              stx.push_back(sim->startPosition().x()); // in disk FF frame
              sty.push_back(sim->startPosition().y());
              stz.push_back(sim->startPosition().z());
           }
           cluSimId_.push_back(sId);
           cluSimPdgId_.push_back(sPdg);
           cluSimCrCode_.push_back(sCr);
           cluSimMom_.push_back(sm);
           cluSimMom2_.push_back(sm2);
           cluSimEdep_.push_back(sed);
           cluSimTime_.push_back(st);
           cluSimPosX_.push_back(spx);
           cluSimPosY_.push_back(spy);
           cluSimPosZ_.push_back(spz);
           cluSimStartX_.push_back(stx);
           cluSimStartY_.push_back(sty);
           cluSimStartZ_.push_back(stz);

           ++nCluster_;
       }


       //--------------------------  Do virtual detectors --------------------------------
       //73/74/77/78 front back inner outer edges disk 0
       //75/76/79/80 front back inner outer edges disk 1
       nVd_ = 0;
       if (vdhits.isValid())
       {
           for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
             {
               const StepPointMC& hit = *iter;

               if ( (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut)
                     && hit.volumeId() != VirtualDetectorId::TT_Back) continue;

                  double hitTime         = fmod(hit.time(),_mbtime);

               CLHEP::Hep3Vector VDPos = cal.geomUtil().mu2eToTracker(hit.position());

               vdId_[nVd_]     = hit.volumeId();
               vdPdgId_[nVd_]  = hit.simParticle()->pdgId();
               vdGenId_[nVd_]  = (hit.simParticle()->genParticle()) ? hit.simParticle()->genParticle()->generatorId().id() : -1;
               vdTime_[nVd_]   = hitTime;
               vdPosX_[nVd_]   = VDPos.x(); //tracker frame
               vdPosY_[nVd_]   = VDPos.y();
               vdPosZ_[nVd_]   = VDPos.z();
               vdMom_[nVd_]    = hit.momentum().mag();
               vdMomX_[nVd_]   = hit.momentum().x();
               vdMomY_[nVd_]   = hit.momentum().y();
               vdMomZ_[nVd_]   = hit.momentum().z();
               vdGenIdx_[nVd_] = hit.simParticle()->generatorIndex();
               ++nVd_;
             }
         }

        Ntup_->Fill();

  }



}

DEFINE_ART_MODULE(mu2e::CaloNNDiag)
