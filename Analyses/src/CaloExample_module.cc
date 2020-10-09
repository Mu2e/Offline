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
#include "fhiclcpp/ParameterSet.h"
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

#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"



namespace mu2e {


  class CaloExample : public art::EDAnalyzer {

     public:

       explicit CaloExample(fhicl::ParameterSet const& pset);
       virtual ~CaloExample() { }

       virtual void beginJob();
       virtual void endJob() {};

       virtual void analyze(const art::Event& e);


     private:

       int diagLevel_;
       bool doGenerated_;

       std::string g4ModuleLabel_;
       std::string generatorModuleLabel_;
       art::InputTag simParticleTag_;


       std::string caloCrystalModuleLabel_;
       std::string caloClusterModuleLabel_;
       std::string caloDigiTruthModuleLabel_;
       std::string caloClusterTruthModuleLabel_;
       std::string virtualDetectorLabel_;
       std::string stepPointMCLabel_;
       SimParticleTimeOffset toff_; 
       int nProcess_;


       TH1F *hcryE_,*hcryT_,*hcryX_,*hcryY_,*hcryZ_;
       TH1F *hcluE_,*hcluT_,*hcluX_,*hcluY_,*hcluZ_,*hcluE_1Et,*hcluE_1E9,*hcluE_1E25,*hcluE_F;       
       TH2F *hxy_,*hECE;

       TTree* Ntup_;
       int   _evt,_run;

       int   nGen_,genPdgId_[16384],genCrCode_[16384];
       float genmomX_[16384],genmomY_[16384],genmomZ_[16384],genStartX_[16384],genStartY_[16384],genStartZ_[16384],genStartT_[16384];

       int   nHits_,cryId_[16384],crySectionId_[16384],crySimIdx_[16384],crySimLen_[16384];
       float cryEtot_,cryTime_[16384],cryEdep_[16384],cryPosX_[16384],cryPosY_[16384],cryPosZ_[16384],_cryLeak[16384];

       int   nSimHit_,crySimId_[16384],crySimPdgId_[16384],crySimCrCode_[16384],crySimGenIdx_[16384],cryConv_[16384];
       float crySimMom_[16384],crySimStartX_[16384],crySimStartY_[16384],crySimStartZ_[16384],crySimStartT_[16384];
       float crySimTime_[16384],crySimEdep_[16348],cryTimeErr_[16348], cryDT_[16384];

       int   nCluster_,nCluSim_,cluNcrys_[16384];
       float cluEnergy_[16384],cluEnergyErr_[16384],cluTime_[16384],cluTimeErr_[16384],cluCogX_[16384],cluCogY_[16384],
             cluCogZ_[16384],cluE1_[16384],cluE9_[16384],cluE25_[16384],cluSecMom_[16384];
       int   cluSplit_[16384],cluConv_[16384],cluSimIdx_[16384],cluSimLen_[16384];
       std::vector<std::vector<int> > cluList_;

       int   cluSimId_[16384],cluSimPdgId_[16384],cluSimGenId_[16384],cluSimGenPdg_[16384],cluSimCrCode_[16384];
       float cluSimMom_[16384],cluSimMom2_[16384],cluSimPosX_[16384],cluSimPosY_[16384],cluSimPosZ_[16384],cluSimStartX_[16384],
             cluSimStartY_[16384],cluSimStartZ_[16384],cluSimTime_[16384],cluSimEdep_[16384];

       int   nVd_,vdId_[16384],vdPdgId_[16384],vdGenId_[16384],vdGenIdx_[16384];
       float vdTime_[16384],vdPosX_[16384],vdPosY_[16384],vdPosZ_[16384],vdMom_[16384],vdMomX_[16384],vdMomY_[16384],vdMomZ_[16384];
  };


//Fix this mess and move to Config
//Add recoDigis

  CaloExample::CaloExample(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    diagLevel_(pset.get<int>("diagLevel",0)),
    doGenerated_(pset.get<bool>("doGenerated",false)),
    g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
    generatorModuleLabel_(pset.get<std::string>("generatorModuleLabel")),
    simParticleTag_(pset.get<std::string>("simParticleTag")),
    caloCrystalModuleLabel_(pset.get<std::string>("caloCrystalModuleLabel")),
    caloClusterModuleLabel_(pset.get<std::string>("caloClusterModuleLabel")),
    caloDigiTruthModuleLabel_(pset.get<std::string>("caloHitTruthModuleLabel")),
    caloClusterTruthModuleLabel_(pset.get<std::string>("caloClusterTruthModuleLabel")),
    virtualDetectorLabel_(pset.get<std::string>("virtualDetectorName")),
    stepPointMCLabel_(pset.get<std::string>("stepPointMCLabel")),
    toff_(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    nProcess_(0),
    Ntup_(0)
  {}

  void CaloExample::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       Ntup_  = tfs->make<TTree>("Calo", "Calo");

       Ntup_->Branch("evt",          &_evt ,         "evt/I");
       Ntup_->Branch("run",          &_run ,         "run/I");
       Ntup_->Branch("cryEtot",      &cryEtot_ ,     "cryEtot/F");

       Ntup_->Branch("nGen",         &nGen_ ,        "nGen/I");
       Ntup_->Branch("genId",        &genPdgId_,     "genId[nGen]/I");
       Ntup_->Branch("genCrCode",    &genCrCode_,    "genCrCode[nGen]/I");
       Ntup_->Branch("genMomX",      &genmomX_,      "genMomX[nGen]/F");
       Ntup_->Branch("genMomY",      &genmomY_,      "genMomY[nGen]/F");
       Ntup_->Branch("genMomZ",      &genmomZ_,      "genMomZ[nGen]/F");
       Ntup_->Branch("genStartX",    &genStartX_,    "genStartX[nGen]/F");
       Ntup_->Branch("genStartY",    &genStartY_,    "genStartY[nGen]/F");
       Ntup_->Branch("genStartZ",    &genStartZ_,    "genStartZ[nGen]/F");
       Ntup_->Branch("genStartT",    &genStartT_,    "genStartT[nGen]/F");

       Ntup_->Branch("nCry",         &nHits_ ,       "nCry/I");
       Ntup_->Branch("cryId",        &cryId_ ,       "cryId[nCry]/I");
       Ntup_->Branch("crySectionId", &crySectionId_, "crySectionId[nCry]/I");
       Ntup_->Branch("cryPosX",      &cryPosX_ ,     "cryPosX[nCry]/F");
       Ntup_->Branch("cryPosY",      &cryPosY_ ,     "cryPosY[nCry]/F");
       Ntup_->Branch("cryPosZ",      &cryPosZ_ ,     "cryPosZ[nCry]/F");
       Ntup_->Branch("cryEdep",      &cryEdep_ ,     "cryEdep[nCry]/F");
       Ntup_->Branch("cryTime",      &cryTime_ ,     "cryTime[nCry]/F");
       Ntup_->Branch("cryTimeErr",   &cryTimeErr_ ,  "cryTimeErr[nCry]/F");
       Ntup_->Branch("cryDT",        &cryDT_ ,       "cryDT[nCry]/F");
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
       Ntup_->Branch("cluCogZ",      &cluCogZ_ ,     "cluCogZ[nCluster]/F");
       Ntup_->Branch("cluNcrys",     &cluNcrys_ ,    "cluNcrys[nCluster]/I");
       Ntup_->Branch("cluE1",        &cluE1_ ,       "cluE1[nCluster]/F");
       Ntup_->Branch("cluE9",        &cluE9_ ,       "cluE9[nCluster]/F");
       Ntup_->Branch("cluE25",       &cluE25_ ,      "cluE25[nCluster]/F");
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


       hcryE_     = tfs->make<TH1F>("cryEdep",  "Energy deposited / crystal", 100,    0., 50.   );
       hcryT_     = tfs->make<TH1F>("cryTime",  "Time of crystal hit",        100,    0., 2000. );
       hcryX_     = tfs->make<TH1F>("cryX",     "X coord of crystal hit",     100,  300., 700.  );
       hcryY_     = tfs->make<TH1F>("cryY",     "Y coord of crystal hit",     100,  300., 700.  );
       hcryZ_     = tfs->make<TH1F>("cryZ",     "Z coord of crystal hit",     100,11000., 13000.);
       hcluE_     = tfs->make<TH1F>("cluEdep",  "Energy deposited / cluster", 150,    0., 150.  );
       hcluT_     = tfs->make<TH1F>("cluTime",  "Time of clustal hit",        100,    0., 2000. );
       hcluX_     = tfs->make<TH1F>("cluX",     "X coord of cluster hit",     100,  300., 700.  );
       hcluY_     = tfs->make<TH1F>("cluY",     "Y coord of cluster hit",     100,  300., 700.  );
       hcluZ_     = tfs->make<TH1F>("cluZ",     "Z coord of cluster hit",     100,11000., 13000.);
       hcluE_1Et  = tfs->make<TH1F>("cluE1Et",  "E1/Etot",                    100,    0., 1.1   );
       hcluE_1E9  = tfs->make<TH1F>("cluE1E9",  "E1/E9",                      100,    0., 1.1   );
       hcluE_1E25 = tfs->make<TH1F>("cluE1E25", "E1/E25",                     100,    0., 1.1   );
       hxy_       = tfs->make<TH2F>("cryxy",    "cry XY",                     350,-700,700,350,-700,700  );

  }



  void CaloExample::analyze(const art::Event& event) 
  {
      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloExample =  "<<nProcess_ <<std::endl;

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      double _mbtime = accPar->deBuncherPeriod;
      toff_.updateMap(event);

      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if (!geom->hasElement<Calorimeter>() ) return;
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(caloCrystalModuleLabel_, caloCrystalHitsHandle);
      const CaloCrystalHitCollection& caloCrystalHits(*caloCrystalHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(caloClusterModuleLabel_, caloClustersHandle);
      const CaloClusterCollection& caloClusters(*caloClustersHandle);

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(g4ModuleLabel_,virtualDetectorLabel_,vdhits);

      //Calo digi truth assignment
      art::Handle<CaloDigiMCTruthAssn> caloDigiTruthHandle;
      event.getByLabel(caloDigiTruthModuleLabel_, caloDigiTruthHandle);
      const CaloDigiMCTruthAssn& caloDigiTruth(*caloDigiTruthHandle);
      
      //Calo cluster truth assignment
      art::Handle<CaloClusterMCTruthAssn> caloClusterTruthHandle;
      event.getByLabel(caloClusterTruthModuleLabel_, caloClusterTruthHandle);
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


       //--------------------------  Do generated particles --------------------------------
       _evt = event.id().event();
       _run = event.run();

       if (diagLevel_ == 3){std::cout << "processing event in calo_example " << nProcess_ << " run and event  = " << _run << " " << _evt << std::endl;}

       
       nGen_=0;
       if (doGenerated_)
       {
           //Get generated particles
           art::Handle<GenParticleCollection> gensHandle;
           event.getByLabel(generatorModuleLabel_, gensHandle);
           const GenParticleCollection& genParticles(*gensHandle);
	   
	   nGen_ = genParticles.size();
	   for (int i=0; i < nGen_; ++i)
	   {
               const GenParticle* gen = &genParticles[i];
               genPdgId_[i]   = gen->pdgId();
               genCrCode_[i]  = gen->generatorId().id();
               genmomX_[i]    = gen->momentum().vect().x();
               genmomY_[i]    = gen->momentum().vect().y();
               genmomZ_[i]    = gen->momentum().vect().z();
               genStartX_[i]  = gen->position().x();
               genStartY_[i]  = gen->position().y();
               genStartZ_[i]  = gen->position().z();
               genStartT_[i]  = gen->time();
	   }
       } 



       //--------------------------  Do calorimeter hits --------------------------------

       nHits_ = nSimHit_ = 0;
       cryEtot_ = 0.0;

       double sumECryConv(0);
       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
       {
           const CaloCrystalHit& hit     = caloCrystalHits.at(ic);
	   int diskId                    = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
            
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
 
           float cryDT(999);
           if (hit.recoCaloDigis().size()==2)
           {
              if (hit.recoCaloDigis().at(0)->ROid()%2==0) cryDT = hit.recoCaloDigis().at(1)->time()- hit.recoCaloDigis().at(0)->time();
              else                                     cryDT = hit.recoCaloDigis().at(0)->time()- hit.recoCaloDigis().at(1)->time();
           } 
 
 
           cryId_[nHits_]        = hit.id();
           crySectionId_[nHits_] = diskId;
           cryEdep_[nHits_]      = hit.energyDep();
           cryTime_[nHits_]      = hit.time();
           cryTimeErr_[nHits_]   = hit.timeErr();
           cryDT_[nHits_]        = cryDT;
           cryPosX_[nHits_]      = crystalPos.x();
           cryPosY_[nHits_]      = crystalPos.y();
           cryPosZ_[nHits_]      = crystalPos.z();
           cryConv_[nHits_]      = isConversion ? 1 : 0;
           cryEtot_             += hit.energyDep();
          
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
            
           hcryE_->Fill(hit.energyDep());
           hcryT_->Fill(hit.time());
           hcryX_->Fill(crystalPos.x());
           hcryY_->Fill(crystalPos.y());
           hcryZ_->Fill(crystalPos.z());
	   hxy_->Fill(crystalPos.x(),crystalPos.y(),hit.energyDep());
           
           if (isConversion) sumECryConv += hit.energyDep();
       }

       //--------------------------  Do clusters --------------------------------
       nCluster_ = nCluSim_ = 0;
       cluList_.clear();
       for (unsigned int ic=0; ic<caloClusters.size();++ic)
       {
          const CaloCluster& cluster = caloClusters.at(ic);
          std::vector<int> cryList;
          for (auto cryPtr : cluster.caloCrystalHitsPtrVector()) cryList.push_back(int(cryPtr.get()- &caloCrystalHits.at(0)));

          auto itMC = caloClusterTruth.begin();
          while (itMC != caloClusterTruth.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}
          unsigned nCluSims   = (itMC != caloClusterTruth.end()) ? itMC->second->nParticles() : 0;

          bool isConversion(false);
          if (nCluSims>0) 
          {
             for (auto& edep : itMC->second->energyDeposits())
             {
                auto parent(edep.sim());
                while (parent->hasParent()) parent = parent->parent();                     
	        if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) isConversion=true;
             }    		          
          }
          
          cluEnergy_[nCluster_]    = cluster.energyDep();
          cluEnergyErr_[nCluster_] = cluster.energyDepErr();
          cluTime_[nCluster_]      = cluster.time();
          cluTimeErr_[nCluster_]   = cluster.timeErr();
          cluNcrys_[nCluster_]     = cluster.size();
          cluCogX_[nCluster_]      = cluster.cog3Vector().x(); //in disk FF frame
          cluCogY_[nCluster_]      = cluster.cog3Vector().y();
          cluCogZ_[nCluster_]      = cluster.cog3Vector().z();
          cluE1_[nCluster_]        = cluster.e1();
          cluE9_[nCluster_]        = cluster.e9();
          cluE25_[nCluster_]       = cluster.e25();
          cluSecMom_[nCluster_]    = cluster.secondMoment();
          cluSplit_[nCluster_]     = cluster.isSplit();
          cluConv_[nCluster_]      = isConversion;
          cluList_.push_back(cryList);

          cluSimIdx_[nCluster_] = nCluSim_;
          cluSimLen_[nCluster_] = nCluSims;

          for (unsigned i=0;i< nCluSims;++i)
	  {	       
	      const auto& eDepMC = itMC->second->energyDeposit(i);	       
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
		 CLHEP::Hep3Vector simPos = cal.geomUtil().mu2eToDiskFF(cluster.diskId(), vdMapEntry->second->position());		  
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

           hcluE_->Fill(cluster.energyDep());
           hcluT_->Fill(cluster.time());
           hcluX_->Fill(cluster.cog3Vector().x());
           hcluY_->Fill(cluster.cog3Vector().y());
           hcluZ_->Fill(cluster.cog3Vector().z());
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
               
               if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;

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
               vdMom_[nVd_]    = hit.momentum().rho();
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

DEFINE_ART_MODULE(mu2e::CaloExample);
