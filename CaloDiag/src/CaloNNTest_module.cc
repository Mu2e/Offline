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

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "CaloCluster/inc/ClusterUtils.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include "Mu2eUtilities/inc/MVATools.hh"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"


namespace 
{
   constexpr int ntupLen = 16384;
}


namespace mu2e {

  class CaloNNTest : public art::EDAnalyzer {

     public:
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Atom<art::InputTag>     caloHitCollection     { Name("caloHitCollection"),      Comment("Calo Hit collection name") }; 
             fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"),  Comment("Calo cluster collection name") }; 
             fhicl::Atom<art::InputTag>     caloHitTruth          { Name("caloHitTruth"),           Comment("CaloHit truth name") }; 
             fhicl::Atom<art::InputTag>     caloClusterTruth      { Name("caloClusterTruth"),       Comment("caloCluster truth name") }; 
             fhicl::Table<MVATools::Config> caloBkgMVA            { Name("caloBkgMVA"),             Comment("MVA Configuration") };
             fhicl::Atom<int>               diagLevel             { Name("diagLevel"),              Comment("Diag Level"),0 };
         };

       explicit CaloNNTest(const art::EDAnalyzer::Table<Config>& config);
       virtual ~CaloNNTest() {}

       virtual void beginJob() override;
       virtual void analyze(const art::Event& e);


     private:
       art::InputTag         caloHitTag_;
       art::InputTag         caloClusterTag_;
       art::InputTag         caloHitTruthTag_;
       art::InputTag         caloClusterTruthTag_;
       MVATools              caloBkgMVA_;
       int                   diagLevel_;
       int                   nProcess_;


       TTree* Ntup_;
       int   _evt,_run;

       int   nHits_,cryId_[ntupLen],crySectionId_[ntupLen],crySimIdx_[ntupLen],crySimLen_[ntupLen];
       float cryTime_[ntupLen],cryEdep_[ntupLen],cryEdepErr_[ntupLen],cryPosX_[ntupLen],cryPosY_[ntupLen],cryPosZ_[ntupLen],_cryLeak[ntupLen];

       int   nSimHit_,crySimId_[ntupLen],crySimPdgId_[ntupLen],crySimCrCode_[ntupLen],crySimGenIdx_[ntupLen],cryConv_[ntupLen];
       float crySimMom_[ntupLen],crySimStartX_[ntupLen],crySimStartY_[ntupLen],crySimStartZ_[ntupLen],crySimStartT_[ntupLen];
       float crySimTime_[ntupLen],crySimEdep_[ntupLen],cryTimeErr_[ntupLen],cryT1_[ntupLen],cryT2_[ntupLen],cryT1Err_[ntupLen],cryT2Err_[ntupLen];

       int   nCluster_,nCluSim_,cluNcrys_[ntupLen];
       float cluEnergy_[ntupLen],cluEnergyErr_[ntupLen],cluTime_[ntupLen],cluTimeErr_[ntupLen],cluCogX_[ntupLen],cluCogY_[ntupLen],cluCogR_[ntupLen],
             cluCogZ_[ntupLen],cluE1_[ntupLen],cluE2_[ntupLen],cluE9_[ntupLen],cluE25_[ntupLen],cluSecMom_[ntupLen],cluEout_[ntupLen],cluEin_[ntupLen];
       int   cluSplit_[ntupLen],cluConv_[ntupLen],cluSimIdx_[ntupLen],cluSimLen_[ntupLen];
       std::vector<std::vector<int> > cluList_;

       int   cluSimId_[ntupLen],cluSimPdgId_[ntupLen],cluSimGenId_[ntupLen],cluSimGenPdg_[ntupLen],cluSimCrCode_[ntupLen];
       float cluSimMom_[ntupLen],cluSimMom2_[ntupLen],cluSimPosX_[ntupLen],cluSimPosY_[ntupLen],cluSimPosZ_[ntupLen],cluSimStartX_[ntupLen],
             cluSimStartY_[ntupLen],cluSimStartZ_[ntupLen],cluSimTime_[ntupLen],cluSimEdep_[ntupLen],cluMVA_[ntupLen];

  };


  CaloNNTest::CaloNNTest(const art::EDAnalyzer::Table<Config>& config) :
    EDAnalyzer{config},
    caloHitTag_         (config().caloHitCollection()),
    caloClusterTag_     (config().caloClusterCollection()),
    caloHitTruthTag_    (config().caloHitTruth()),
    caloClusterTruthTag_(config().caloClusterTruth()),
    caloBkgMVA_         (config().caloBkgMVA()),
    diagLevel_          (config().diagLevel()),
    nProcess_(0),
    Ntup_(0)
  {}

  void CaloNNTest::beginJob(){

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
       Ntup_->Branch("cluEout",      &cluEout_ ,     "cluEout[nCluster]/F");
       Ntup_->Branch("cluEin",       &cluEin_ ,      "cluEin[nCluster]/F");
       Ntup_->Branch("cluSecMom",    &cluSecMom_ ,   "cluSecMom[nCluster]/F");
       Ntup_->Branch("cluSplit",     &cluSplit_ ,    "cluSplit[nCluster]/I");
       Ntup_->Branch("cluMVA",       &cluMVA_ ,      "cluMVA[nCluster]/F");
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
       Ntup_->Branch("cluSimStartX", &cluSimStartX_ ,"cluSimStartX[nCluSim]/F");
       Ntup_->Branch("cluSimStartY", &cluSimStartY_ ,"cluSimStartY[nCluSim]/F");
       Ntup_->Branch("cluSimStartZ", &cluSimStartZ_ ,"cluSimStartZ[nCluSim]/F");
       Ntup_->Branch("cluSimTime",   &cluSimTime_ ,  "cluSimTime[nCluSim]/F");
       Ntup_->Branch("cluSimEdep",   &cluSimEdep_ ,  "cluSimEdep[nCluSim]/F");

       caloBkgMVA_.initMVA();
  }

  void CaloNNTest::analyze(const art::Event& event) 
  {
      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloNNTest =  "<<nProcess_ <<std::endl;


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

      //Calo digi truth assignment
      art::Handle<CaloHitMCCollection> caloHitMCHandle;
      event.getByLabel(caloHitTruthTag_, caloHitMCHandle);
      const CaloHitMCCollection& caloHitTruth(*caloHitMCHandle);
      
      //Calo cluster truth assignment
      art::Handle<CaloClusterMCCollection> caloClusterMCHandle;
      event.getByLabel(caloClusterTruthTag_, caloClusterMCHandle);
      const CaloClusterMCCollection& caloClusterTruth(*caloClusterMCHandle);

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

         const auto eDepMCs =  caloHitTruth[ic].energyDeposits();

         bool isConversion(false); 
         for (auto& edep : eDepMCs)
         {
            if (edep.sim()->creationCode() == ProcessCode::mu2eCeMinusEndpoint)  isConversion=true;
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
         crySimLen_[nHits_]    = eDepMCs.size();

         double sumEdepMC(0),edepTime(0);
         for (unsigned i=0;i< eDepMCs.size();++i)
	 {	       
	     const auto& eDepMC = eDepMCs[i];	       

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

             sumEdepMC += eDepMC.energyDep();
             if (edepTime<1) edepTime = eDepMC.time();
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
           const auto& cluster = caloClusters.at(ic);
           const auto& eDepMCs = caloClusterTruth[ic].energyDeposits();

           for (auto& edep : eDepMCs)
           {
               if (edep.sim()->creationCode() == ProcessCode::mu2eCeMinusEndpoint && cluster.energyDep() > convEnergy) 
               {
		  convEnergy = cluster.energyDep();
		  icMCIdx = ic;
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

          const auto& eDepMCs =  caloClusterTruth[ic].energyDeposits();
          //bool isConversion(false); 
          //for (auto& edep : eDepMCs)
          //{
          //   if (edep.sim()->creationCode() == ProcessCode::mu2eCeMinusEndpoint)  isConversion=true;
          //}    		          

          
          //If cliuster energy is below some threshold (e.g. 40 MeV), then simply put the MVA score to zero and skip this part
	  const CaloHit& seedHit    = CaloHits.at(cryList[0]);
          CLHEP::Hep3Vector seedPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(),cal.crystal(seedHit.crystalID()).position());
	  const auto& neighborsId   = cal.crystal(seedHit.crystalID()).neighbors();
	  const auto& nneighborsId  = cal.crystal(seedHit.crystalID()).nextNeighbors();
	  

          double enerIn(0),enerOut(0),e2(seedHit.energyDep()), e9(seedHit.energyDep()),e25(seedHit.energyDep());
          double r0 = seedPos.perp();
          for (auto cryPtr : cluster.caloHitsPtrVector())
	  {
              const CaloHit& hit    = *cryPtr;
              CLHEP::Hep3Vector pos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(),cal.crystal(hit.crystalID()).position());
	      double r1 = pos.perp();
	      if (r1 > 1.01*r0) enerOut += hit.energyDep();	     
	      if (r1 < 0.99*r0) enerIn  += hit.energyDep();	     
              
	      if (std::find(neighborsId.begin(),  neighborsId.end(),  hit.crystalID()) != neighborsId.end())  {e9 += hit.energyDep();e25 += hit.energyDep();}
	      if (std::find(nneighborsId.begin(), nneighborsId.end(), hit.crystalID()) != nneighborsId.end()) {e25 += hit.energyDep();}
	  }
	  
	  if (cryList.size()>1) e2 += CaloHits.at(cryList[1]).energyDep();


	  std::vector<float> mvavars(8,0.0);
	  mvavars[0] = cluster.energyDep();
	  mvavars[1] = cluster.cog3Vector().perp();
	  mvavars[2] = cluster.size(); 
	  mvavars[3] = seedHit.energyDep();
	  mvavars[4] = e2;
	  mvavars[5] = e9;
	  mvavars[6] = e25;
	  mvavars[7] = cluster.diskID();
	  
	  float mvaout = caloBkgMVA_.evalMVA(mvavars);


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
          cluEout_[nCluster_]      = enerOut;
          cluEin_[nCluster_]       = enerIn;
          cluSecMom_[nCluster_]    = cluUtil.secondMoment();
          cluSplit_[nCluster_]     = cluster.isSplit();
          cluMVA_[nCluster_]       = mvaout;
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

              cluSimId_[nCluSim_]     = sim->id().asInt();
              cluSimPdgId_[nCluSim_]  = sim->pdgId();
              cluSimCrCode_[nCluSim_] = sim->creationCode();
              cluSimGenId_[nCluSim_]  = genId;
	      cluSimGenPdg_[nCluSim_] = genPdg;
              cluSimTime_[nCluSim_]   = eDepMC.time();
              cluSimEdep_[nCluSim_]   = eDepMC.energyDep();
              cluSimMom_[nCluSim_]    = eDepMC.momentumIn();
              cluSimStartX_[nCluSim_] = sim->startPosition().x(); // in disk FF frame
              cluSimStartY_[nCluSim_] = sim->startPosition().y();
              cluSimStartZ_[nCluSim_] = sim->startPosition().z();  

              ++nCluSim_;
           }
           ++nCluster_;
       }

       Ntup_->Fill();

  }



}  

DEFINE_ART_MODULE(mu2e::CaloNNTest);
