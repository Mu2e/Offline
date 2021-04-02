//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
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

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "CaloCluster/inc/ClusterUtils.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

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
             using SPTO    = SimParticleTimeOffset::Config;
             fhicl::Atom<art::InputTag>     vdCollection          { Name("vdCollection"),           Comment("Virtual detector collection name") }; 
             fhicl::Table<SPTO>             timeOffsets           { Name("timeOffsets"),            Comment("Sim Particle Time Offset Maps")};
             fhicl::Atom<art::InputTag>     caloHitCollection     { Name("caloHitCollection"),      Comment("Calo Hit collection name") }; 
             fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"),  Comment("Calo cluster collection name") }; 
             fhicl::Atom<art::InputTag>     caloHitTruth          { Name("caloHitTruth"),           Comment("CaloHit truth name") }; 
             fhicl::Atom<art::InputTag>     caloClusterTruth      { Name("caloClusterTruth"),       Comment("caloCluster truth name") }; 
             fhicl::Atom<int>               diagLevel             { Name("diagLevel"),              Comment("Diag Level"),0 };
         };

       explicit CaloNNDiag(const art::EDAnalyzer::Table<Config>& config);
       virtual ~CaloNNDiag() {}

       virtual void beginJob();
       virtual void endJob() {};
       virtual void analyze(const art::Event& e);


     private:
       art::InputTag         virtualDetectorTag_;
       SimParticleTimeOffset toff_; 
       art::InputTag         caloHitTag_;
       art::InputTag         caloClusterTag_;
       art::InputTag         caloHitTruthTag_;
       art::InputTag         caloClusterTruthTag_;
       int                   diagLevel_;
       int                   nProcess_;


       TTree* Ntup_;
       int   _evt,_run;

       int   nHits_,cryId_[ntupLen],crySectionId_[ntupLen],crySimIdx_[ntupLen],crySimLen_[ntupLen];
       float cryTime_[ntupLen],cryEdep_[ntupLen],cryEdepErr_[ntupLen],cryPosX_[ntupLen],cryPosY_[ntupLen],cryPosZ_[ntupLen],_cryLeak[ntupLen];

       int   nSimHit_,crySimId_[ntupLen],crySimPdgId_[ntupLen],crySimCrCode_[ntupLen],crySimGenIdx_[ntupLen],cryConv_[ntupLen];
       float crySimMom_[ntupLen],crySimStartX_[ntupLen],crySimStartY_[ntupLen],crySimStartZ_[ntupLen],crySimStartT_[ntupLen];
       float crySimTime_[ntupLen],crySimEdep_[ntupLen],cryTimeErr_[ntupLen],cryT1_[ntupLen],cryT2_[ntupLen],cryT1Err_[ntupLen],cryT2Err_[ntupLen];

       int   nCluster_,nCluSim_,cluNcrys_[ntupLen],cluDisk_[ntupLen];
       float cluEnergy_[ntupLen],cluEnergyErr_[ntupLen],cluTime_[ntupLen],cluTimeErr_[ntupLen],cluCogX_[ntupLen],cluCogY_[ntupLen],cluCogR_[ntupLen],
             cluCogZ_[ntupLen],cluE1_[ntupLen],cluE2_[ntupLen],cluE9_[ntupLen],cluE25_[ntupLen],cluSecMom_[ntupLen],cluDR_[ntupLen];
       int   cluSplit_[ntupLen],cluConv_[ntupLen],cluSimIdx_[ntupLen],cluSimLen_[ntupLen];
       std::vector<std::vector<int> > cluList_;

       int   cluSimId_[ntupLen],cluSimPdgId_[ntupLen],cluSimGenId_[ntupLen],cluSimGenPdg_[ntupLen],cluSimCrCode_[ntupLen];
       float cluSimMom_[ntupLen],cluSimMom2_[ntupLen],cluSimPosX_[ntupLen],cluSimPosY_[ntupLen],cluSimPosZ_[ntupLen],cluSimStartX_[ntupLen],
             cluSimStartY_[ntupLen],cluSimStartZ_[ntupLen],cluSimTime_[ntupLen],cluSimEdep_[ntupLen];

       int   nVd_,vdId_[ntupLen],vdPdgId_[ntupLen],vdGenId_[ntupLen],vdGenIdx_[ntupLen];
       float vdTime_[ntupLen],vdPosX_[ntupLen],vdPosY_[ntupLen],vdPosZ_[ntupLen],vdMom_[ntupLen],vdMomX_[ntupLen],vdMomY_[ntupLen],vdMomZ_[ntupLen];
  };


  CaloNNDiag::CaloNNDiag(const art::EDAnalyzer::Table<Config>& config) :
    EDAnalyzer{config},
    virtualDetectorTag_ (config().vdCollection()),
    toff_               (config().timeOffsets()),
    caloHitTag_         (config().caloHitCollection()),
    caloClusterTag_     (config().caloClusterCollection()),
    caloHitTruthTag_    (config().caloHitTruth()),
    caloClusterTruthTag_(config().caloClusterTruth()),
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

       Ntup_->Branch("crySimIdx",    &crySimIdx_ ,   "crySimIdx[nCry]/I");
       Ntup_->Branch("crySimLen",    &crySimLen_ ,   "crySimLen[nCry]/I");
       Ntup_->Branch("nSim",         &nSimHit_ ,     "nSim/I");
       Ntup_->Branch("simId",        &crySimId_ ,    "simId[nSim]/I");
       Ntup_->Branch("simPdgId",     &crySimPdgId_ , "simPdgId[nSim]/I");
       Ntup_->Branch("simCrCode",    &crySimCrCode_ ,"simCrCode[nSim]/I");
       Ntup_->Branch("simMom",       &crySimMom_ ,   "simMom[nSim]/F");
       Ntup_->Branch("simStartX",    &crySimStartX_ ,"simStartX[nSim]/F");
       Ntup_->Branch("simStartY",    &crySimStartY_ ,"simStartY[nSim]/F");
       Ntup_->Branch("simStartZ",    &crySimStartZ_ ,"simStartZ[nSim]/F");
       Ntup_->Branch("simStartT",    &crySimStartT_ ,"simStartT[nSim]/F");
       Ntup_->Branch("simTime",      &crySimTime_ ,  "simTime[nSim]/F");
       Ntup_->Branch("simEdep",      &crySimEdep_ ,  "simEdep[nSim]/F");
       Ntup_->Branch("simGenIdx",    &crySimGenIdx_ ,"simGenIdx[nSim]/I");

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
       Ntup_->Branch("cluDR",        &cluDR_ ,       "cluDR[nCluster]/F");
       Ntup_->Branch("cluSecMom",    &cluSecMom_ ,   "cluSecMom[nCluster]/F");
       Ntup_->Branch("cluSplit",     &cluSplit_ ,    "cluSplit[nCluster]/I");
       Ntup_->Branch("cluConv",      &cluConv_ ,     "cluConv[nCluster]/I");
       Ntup_->Branch("cluSimIdx",    &cluSimIdx_ ,   "cluSimIdx[nCluster]/I");
       Ntup_->Branch("cluSimLen",    &cluSimLen_ ,   "cluSimLen[nCluster]/I");
       Ntup_->Branch("cluList",      &cluList_);

       Ntup_->Branch("nCluSim",      &nCluSim_ ,     "nCluSim/I");
       Ntup_->Branch("cluSimId",     &cluSimId_ ,    "cluSimId[nCluSim]/I");
       Ntup_->Branch("cluSimPdgId",  &cluSimPdgId_ , "cluSimPdgId[nCluSim]/I");
       Ntup_->Branch("cluSimGenId",  &cluSimGenId_ , "cluSimGenId[nCluSim]/I");
       Ntup_->Branch("cluSimGenPdg", &cluSimGenPdg_, "cluSimGenPdg[nCluSim]/I");
       Ntup_->Branch("cluSimCrCode", &cluSimCrCode_ ,"cluSimCrCode[nCluSim]/I");
       Ntup_->Branch("cluSimMom",    &cluSimMom_ ,   "cluSimMom[nCluSim]/F");
       Ntup_->Branch("cluSimMom2",   &cluSimMom2_ ,  "cluSimMom2[nCluSim]/F");
       Ntup_->Branch("cluSimPosX",   &cluSimPosX_ ,  "cluSimPosX[nCluSim]/F");
       Ntup_->Branch("cluSimPosY",   &cluSimPosY_ ,  "cluSimPosY[nCluSim]/F");
       Ntup_->Branch("cluSimPosZ",   &cluSimPosZ_ ,  "cluSimPosZ[nCluSim]/F");
       Ntup_->Branch("cluSimStartX", &cluSimStartX_ ,"cluSimStartX[nCluSim]/F");
       Ntup_->Branch("cluSimStartY", &cluSimStartY_ ,"cluSimStartY[nCluSim]/F");
       Ntup_->Branch("cluSimStartZ", &cluSimStartZ_ ,"cluSimStartZ[nCluSim]/F");
       Ntup_->Branch("cluSimTime",   &cluSimTime_ ,  "cluSimTime[nCluSim]/F");
       Ntup_->Branch("cluSimEdep",   &cluSimEdep_ ,  "cluSimEdep[nCluSim]/F");

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

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      double _mbtime = accPar->deBuncherPeriod;
      toff_.updateMap(event);

      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if (!geom->hasElement<Calorimeter>() ) return;
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(caloHitTag_, CaloHitsHandle);
      const CaloHitCollection& CaloHits(*CaloHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(caloClusterTag_, caloClustersHandle);
      const CaloClusterCollection& caloClusters(*caloClustersHandle);

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(virtualDetectorTag_,vdhits);

      //Calo digi truth assignment
      art::Handle<CaloHitMCTruthAssn> caloDigiTruthHandle;
      event.getByLabel(caloHitTruthTag_, caloDigiTruthHandle);
      const CaloHitMCTruthAssn& caloDigiTruth(*caloDigiTruthHandle);
      
      //Calo cluster truth assignment
      art::Handle<CaloClusterMCTruthAssn> caloClusterTruthHandle;
      event.getByLabel(caloClusterTruthTag_, caloClusterTruthHandle);
      const CaloClusterMCTruthAssn& caloClusterTruth(*caloClusterTruthHandle);




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

       nHits_ = nSimHit_ = 0;

       for (unsigned int ic=0; ic<CaloHits.size();++ic)
       {
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
              int idx0 = cal.caloIDMapper().SiPMIdx(hit.recoCaloDigis().at(0)->SiPMID());
              int idx1 = cal.caloIDMapper().SiPMIdx(hit.recoCaloDigis().at(1)->SiPMID());
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
          
           crySimIdx_[nHits_]    = nSimHit_;
           crySimLen_[nHits_]    = nCrySims;
           
           for (unsigned i=0;i< nCrySims;++i)
	   {	                      
	       const auto& eDepMC = itMC->second->energyDeposit(i);
               
               auto parent(eDepMC.sim());
               while (parent->hasParent()) parent = parent->parent();               
               int genId=-1;
               if (parent->genParticle()) genId = parent->genParticle()->generatorId().id();
              
	       crySimId_[nSimHit_]      = eDepMC.sim()->id().asInt();
               crySimPdgId_[nSimHit_]   = eDepMC.sim()->pdgId();
               crySimCrCode_[nSimHit_]  = eDepMC.sim()->creationCode();
       	       crySimTime_[nSimHit_]    = eDepMC.time();
               crySimEdep_[nSimHit_]    = eDepMC.energyDep();	                      
               crySimMom_[nSimHit_]     = eDepMC.momentumIn();	       
	       crySimStartX_[nSimHit_]  = parent->startPosition().x();
	       crySimStartY_[nSimHit_]  = parent->startPosition().y();
	       crySimStartZ_[nSimHit_]  = parent->startPosition().z();
	       crySimStartT_[nSimHit_]  = parent->startGlobalTime();
	       crySimGenIdx_[nSimHit_]  = genId;
	       ++nSimHit_;
           }
           ++nHits_;            
       }

       //--------------------------  Do clusters --------------------------------
       nCluster_ = nCluSim_ = 0;
       cluList_.clear();


       //find the most energetic CE cluster
       unsigned icMCIdx(-1); double convEnergy(0);
       for (unsigned  ic=0; ic<caloClusters.size();++ic)
       {
           const CaloCluster& cluster = caloClusters.at(ic);
           auto itMC = caloClusterTruth.begin();
           while (itMC != caloClusterTruth.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}

           if (itMC != caloClusterTruth.end()) 
           {
               for (auto& edep : itMC->second->energyDeposits())
               {
                   auto parent(edep.sim());
                   while (parent->hasParent()) parent = parent->parent();                     
	           if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() && cluster.energyDep() > convEnergy)
		   {
		      convEnergy = cluster.energyDep();
		      icMCIdx = ic;
		   } 
               }    		          
           }
       }



       for (unsigned  ic=0; ic<caloClusters.size();++ic)
       {
          const CaloCluster& cluster = caloClusters.at(ic);
          std::vector<int> cryList;
          for (auto cryPtr : cluster.caloHitsPtrVector()) cryList.push_back(std::distance(&CaloHits.at(0),cryPtr.get()));

          ClusterUtils cluUtil(cal, cluster);
          auto cog = cluUtil.cog3Vector();

          auto itMC = caloClusterTruth.begin();
          while (itMC != caloClusterTruth.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}
          const auto eDepMCs = (itMC != caloClusterTruth.end()) ? itMC->second->energyDeposits() : std::vector<CaloEDepMC>{};
          
	  const CaloHit& seedHit    = CaloHits.at(cryList[0]);
          CLHEP::Hep3Vector seedPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(),cal.crystal(seedHit.crystalID()).position());

          double dr(0);
          if (cryList.size()>1) 
	  {
              const CaloHit& hit    = CaloHits.at(cryList[1]);
              CLHEP::Hep3Vector hitPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(),cal.crystal(hit.crystalID()).position());
	      double dx = hitPos.x()-seedPos.x();
              double dy = hitPos.y()-seedPos.y();
              dr = sqrt(dx*dx+dy*dy);	     
	  }


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
          cluDR_[nCluster_]        = dr;
          cluSecMom_[nCluster_]    = cluUtil.secondMoment();
          cluSplit_[nCluster_]     = cluster.isSplit();
          cluConv_[nCluster_]      = ic == icMCIdx;
          cluList_.push_back(cryList);

          cluSimIdx_[nCluster_] = nCluSim_;
          cluSimLen_[nCluster_] = eDepMCs.size();

          for (unsigned i=0;i< eDepMCs.size();++i)
	  {	       
	      const auto& eDepMC = eDepMCs[i];	       
              art::Ptr<SimParticle> sim = eDepMC.sim();

	      art::Ptr<SimParticle> smother(sim);
              while (smother->hasParent() && !smother->genParticle() ) smother = smother->parent();
              int genId=-1;
              if (smother->genParticle()) genId = smother->genParticle()->generatorId().id();
              int genPdg=-1;
              if (smother->genParticle()) genPdg = smother->genParticle()->pdgId();


	      double simMom(-1);
	      CLHEP::Hep3Vector simPos(0,0,0);
	      auto vdMapEntry = vdMap.find(sim);
	      if (vdMapEntry != vdMap.end())
	      {
	         simMom = vdMapEntry->second->momentum().mag();
		 CLHEP::Hep3Vector simPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(), vdMapEntry->second->position());		  
	      } 

              cluSimId_[nCluSim_]     = sim->id().asInt();
              cluSimPdgId_[nCluSim_]  = sim->pdgId();
              cluSimCrCode_[nCluSim_] = sim->creationCode();
              cluSimGenId_[nCluSim_]  = genId;
	      cluSimGenPdg_[nCluSim_] = genPdg;
              cluSimTime_[nCluSim_]   = eDepMC.time();
              cluSimEdep_[nCluSim_]   = eDepMC.energyDep();
              cluSimMom_[nCluSim_]    = eDepMC.momentumIn();
              cluSimMom2_[nCluSim_]   = simMom;
              cluSimPosX_[nCluSim_]   = simPos.x(); // in disk FF frame
              cluSimPosY_[nCluSim_]   = simPos.y();
              cluSimPosZ_[nCluSim_]   = simPos.z();  
              cluSimStartX_[nCluSim_] = sim->startPosition().x(); // in disk FF frame
              cluSimStartY_[nCluSim_] = sim->startPosition().y();
              cluSimStartZ_[nCluSim_] = sim->startPosition().z();  

              ++nCluSim_;
           }
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
  
               double hitTimeUnfolded = toff_.totalTimeOffset(hit.simParticle()) + hit.time();
   	       double hitTime         = fmod(hitTimeUnfolded,_mbtime);

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

DEFINE_ART_MODULE(mu2e::CaloNNDiag);
