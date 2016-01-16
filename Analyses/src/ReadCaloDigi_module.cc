//
// An EDAnalyzer module that reads back the hits created by the Calorimeter Digitization chain
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author 

// ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TSpline.h"
#include "TFile.h"
#include "CalPatRec/inc/THackData.hh"
#include "TROOT.h"
#include "TFolder.h"
#include "TTree.h"
#include "TH2F.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"


#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/RecoCaloDigi.hh"
#include "RecoDataProducts/inc/RecoCaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>



using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;


namespace mu2e {
  
  class ReadCaloDigi : public art::EDAnalyzer {
     
  public:

    explicit ReadCaloDigi(fhicl::ParameterSet const& pset);
    virtual ~ReadCaloDigi() { }

    virtual void beginRun(art::Run& );

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);

       
  private:
       
    typedef std::vector< art::Handle<mu2e::StepPointMCCollection> > StepMCHandleVector;
    typedef art::Ptr<mu2e::CaloCrystalHit>                          CaloCrystalHitPtr;

    void             initVHits     ();

    int                        _diagLevel; 
    
    std::string                _caloDigisModuleLabel;
    std::string                _caloCrystalModuleLabel;
    std::string                _stepPoints;
    std::string                _rostepPoints;
    std::string                _caloClusterModuleLabel;
    std::string                _vdStepPoints;

    SimParticleTimeOffset      _toff;     // time offset smearing
    double                     _mbtime;
    double                     _mbbuffer;
    double                     _blindTime;

    double                     _psdThreshold;

    TTree*                     _Ntup;
    int                        _nProcess;

    int                        _evt,_run,_caloCrystals,_caloDisk0Crystals,_caloDisk1Crystals,_nHits,_nCluster, _nCryRO;
       
    int                        _cryId[16384],_crySectionId[16384];

    int                        _cluNcrys[204828];
    
    float                      _caloVolume, _crystalVolume;

    float                      _cryEtot,_cryTime[16384],_cryEdep[16384],_cryDose[16384];
    float                      _cryMCTime    [16384];
    float                      _cryMCEdep    [16384];
    int                        _cryNParticles[16384];
    float                      _cryPsd       [16384];
    float                      _cryIsConv    [16384];

    float                      _cryPosX[16384],_cryPosY[16384],_cryPosZ[16384], _cryPosR[16384];
    
    float                      _cluEnergy[2048], _cluCrysE[2048], _cluTime[2048], _cluMCTime[2048], _cluMCMeanTime[2048];
    float                      _cluCogX[2048],_cluCogY[2048],_cluCogZ[2048], _cluCogR[2048];

    int                        _cluConv[2048];

    Int_t                     _vNHits;
    Float_t                   _vP, _vPx, _vPy, _vPz, _vPt;
    Float_t                   _vE, _vEKin, _vM;
    Float_t                   _vT;
    Float_t                   _vX, _vY, _vZ;
    Float_t                   _vCosth, _vRadius;
    Int_t                     _vId;
    Int_t                     _vPdgId;

    const Calorimeter*        _calorimeter; // cached pointer to the calorimeter geometry

  };


  void     ReadCaloDigi::initVHits(){

    _vNHits  = 0;
      
    _vP      = 0.;
    _vPx     = 0.;
    _vPy     = 0.;
    _vPz     = 0.;
    _vPt     = 0.;
    _vPdgId  = 0.;
    _vM      = 0.;
    _vE      = 0.;
    _vEKin   = 0.;
    _vT      = 0.;
    _vX      = 0.;
    _vY      = 0.;
    _vZ      = 0.;
          
    _vCosth  = 0.;
    _vRadius = 0.;
    _vId     = 0.;
  }
  
  
  ReadCaloDigi::ReadCaloDigi(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel                     (pset.get<int>   ("diagLevel")),	     
    _caloDigisModuleLabel          (pset.get<string>("caloDigisModuleLabel")),
    _caloCrystalModuleLabel        (pset.get<string>("caloCrystalModuleLabel")),
    _stepPoints                    (pset.get<string>("calorimeterStepPoints")),
    _rostepPoints                  (pset.get<string>("calorimeterROStepPoints")),
    _caloClusterModuleLabel        (pset.get<string>("caloClusterModuleLabel")),
    _vdStepPoints                  (pset.get<string>("vdStepPoints")),	     
    _toff                          (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _mbbuffer                      (pset.get<double>             ("TimeFoldingBuffer")),  // ns
    _blindTime                     (pset.get<double>             ("blindTime" )),         // ns
    _psdThreshold                  (pset.get<double>("psdThreshold")),        
    _Ntup(0),_nProcess(0)

  {}

  void ReadCaloDigi::beginRun(art::Run& ){
//     mu2e::GeomHandle<mu2e::Calorimeter> ch;
//     _calorimeter = ch.get();
  }
  void ReadCaloDigi::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("Calo", "Calo");



    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("caloCrystals", &_caloCrystals ,"caloCrystals/I");
    _Ntup->Branch("caloDisk0Crystals", &_caloDisk0Crystals ,"caloDisk0Crystals/I");
    _Ntup->Branch("caloDisk1Crystals", &_caloDisk1Crystals ,"caloDisk1Crystals/I");
     
    _Ntup->Branch("cryEtot",      &_cryEtot ,    "cryEtot/F");

    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryPosR",      &_cryPosR ,     "cryPosR[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryMCTime",    &_cryMCTime    ,"cryMCTime[nCry]/F");
    _Ntup->Branch("cryMCEdep",    &_cryMCEdep    ,"cryMCEdep[nCry]/F");
    _Ntup->Branch("cryNParticles",&_cryNParticles,"cryNParticles[nCry]/I");
    _Ntup->Branch("cryPsd",       &_cryPsd       ,"cryPsd[nCry]/F");
    _Ntup->Branch("cryIsConv",    &_cryIsConv    ,"cryIsConv[nCry]/I");
    _Ntup->Branch("nCluster",     &_nCluster ,    "nCluster/I");
    _Ntup->Branch("cluEnergy",    &_cluEnergy ,   "cluEnergy[nCluster]/F");
    _Ntup->Branch("cluCrysE",     &_cluCrysE ,    "cluCrysE[nCluster]/F");
    _Ntup->Branch("cluTime",      &_cluTime ,     "cluTime[nCluster]/F");
    _Ntup->Branch("cluMCTime",    &_cluMCTime ,   "cluMCTime[nCluster]/F");
    _Ntup->Branch("cluMCMeanTime",    &_cluMCMeanTime ,   "cluMCMeanTime[nCluster]/F");
    _Ntup->Branch("cluCogX",      &_cluCogX ,     "cluCogX[nCluster]/F");	
    _Ntup->Branch("cluCogY",      &_cluCogY ,     "cluCogY[nCluster]/F");	
    _Ntup->Branch("cluCogZ",      &_cluCogZ ,     "cluCogZ[nCluster]/F");	
    _Ntup->Branch("cluCogR",      &_cluCogR ,     "cluCogR[nCluster]/F");	
    _Ntup->Branch("cluNcrys",     &_cluNcrys ,    "cluNCrys[nCluster]/I");	
    _Ntup->Branch("cluConv",      &_cluConv ,     "cluConv[nCluster]/I");	
  
    _Ntup->Branch("vNHits",   &_vNHits ,  "vNHits/I");
    _Ntup->Branch("vId"   ,   &_vId ,  "vId/I");
    _Ntup->Branch("vPdgId",   &_vPdgId ,  "vPdgId/I");

    _Ntup->Branch("vP" ,   &_vP ,  "vP/F");
    _Ntup->Branch("vPx",   &_vPx ,  "vPx/F");
    _Ntup->Branch("vPy",   &_vPy ,  "vPy/F");
    _Ntup->Branch("vPz",   &_vPz ,  "vPz/F");
    _Ntup->Branch("vE" ,   &_vE ,  "vE/F");
    _Ntup->Branch("vEKin",   &_vEKin ,  "vEKin/F");
    _Ntup->Branch("vM",   &_vM ,  "vM/F");
    _Ntup->Branch("vT",   &_vT ,  "vT/F");
    _Ntup->Branch("vX",   &_vX ,  "vX/F");
    _Ntup->Branch("vY",   &_vY ,  "vY/F");
    _Ntup->Branch("vZ",   &_vZ ,  "vZ/F");
    _Ntup->Branch("vCosth",   &_vCosth ,  "vCosth/F");
    _Ntup->Branch("vRadius",   &_vRadius ,  "vRadius/F");



  }



  void ReadCaloDigi::endJob(){
  }




  void ReadCaloDigi::analyze(const art::Event& event) {
    if (_nProcess == 0){
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      _calorimeter = ch.get();
    }

    ++_nProcess;

    //load the timeoffset
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);

    //data about hits in the calorimeter crystals
    art::Handle<CaloDigiMCCollection> caloDigiMCHandle;
    event.getByLabel(_caloDigisModuleLabel, caloDigiMCHandle);
    if (!caloDigiMCHandle.isValid()){
      printf("[ReadCaloDigi::analyze] no CaloDigiMCCollection: BAILS OUT\n");
      return;
    }
    const     CaloDigiMCCollection*       caloDigiMCCol = caloDigiMCHandle.product();
    

    art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
    event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
    if (!caloCrystalHitsHandle.isValid()){
      printf("[ReadCaloDigi::analyze] no CaloCrystalHitCollection: BAILS OUT \n");
    }
    const     CaloCrystalHitCollection*  caloCrystalHits = caloCrystalHitsHandle.product();


    //data about clusters
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    if (!caloClustersHandle.isValid()){
      printf("[ReadCaloDigi::analyze] no CaloClusterCollection: BAILS OUT \n");
    }
    const     CaloClusterCollection*  caloClusters = caloClustersHandle.product();


    //Handle to VD steps
    art::ProductInstanceNameSelector selector_vdhits("virtualdetector");
    StepMCHandleVector vdStepsHandleVec;
    art::Handle<StepPointMCCollection> *vdStepsHandle;
    const StepPointMCCollection *vdHits;
    event.getMany(selector_vdhits, vdStepsHandleVec);


    GlobalConstantsHandle<ParticleDataTable> pdt;

  
    _evt          = event.id().event();
    _run          = event.run();
    _caloCrystals = _calorimeter->nCrystal();
    _caloDisk0Crystals = _calorimeter->section(0).nCrystals();
    _caloDisk1Crystals = _calorimeter->section(1).nCrystals();
    
    int vdVecSize = vdStepsHandleVec.size();
    
    initVHits();

    for (int j=0; j< vdVecSize; ++j){
      vdStepsHandle = & vdStepsHandleVec[j];
      if (vdStepsHandle->isValid())  {
	vdHits        = vdStepsHandle->operator->();

	for (size_t i=0; i<vdHits->size(); ++i) {
	  StepPointMC hit = vdHits->at(i);

	  if (hit.simParticle()->fromGenerator()) {
	    int id = hit.volumeId();
	    
	    if (id == VirtualDetectorId::EMC_Disk_0_SurfIn  ||
		id == VirtualDetectorId::EMC_Disk_1_SurfIn  ||
		id == VirtualDetectorId::EMC_Disk_0_EdgeIn  ||
		id == VirtualDetectorId::EMC_Disk_1_EdgeIn    ) {


	      art::Ptr<SimParticle> const& simptr = hit.simParticle();
	      //2016-01-10 G. PEzzullo temporary comment for using 
	      // a custom made gen particle

	      SimParticle const& sim  = *simptr;
	      if ( sim.fromGenerator() ){
	      	GenParticle* gen = (GenParticle*) &(*sim.genParticle());
	      	if ( gen->generatorId() != GenId::conversionGun ){
	      	  continue;
	      	}
	      }

	      _vP      = hit.momentum().mag();
	      _vPx     = hit.momentum().x();
	      _vPy     = hit.momentum().y();
	      _vPz     = hit.momentum().z();
	      _vPt     = std::sqrt( std::pow(_vPx,2.)+std::pow(_vPy,2.) );
	      _vPdgId  = hit.simParticle()->pdgId();
	      _vM      = pdt->particle(_vPdgId).ref().mass();
	      _vE      = sqrt(_vP*_vP + _vM*_vM);
	      _vEKin   = _vE - _vM;
	      

	      double     hitTime = fmod(hit.time() + _toff.totalTimeOffset(simptr), _mbtime);
	      if (hitTime < _mbbuffer) {
		if (hitTime+_mbtime > _blindTime) {
		  hitTime = hitTime + _mbtime;
		}
	      }
	      else {
		if (hitTime > (_mbtime - _mbbuffer)) {
		  if (hitTime - _mbtime > _blindTime) {
		    hitTime =   hitTime - _mbtime;
		  }
		}
	      }

	      _vT      = hitTime;
	      _vX      = hit.position().x()+3904.;
	      _vY      = hit.position().y();
	      _vZ      = hit.position().z();
	                
	      _vCosth  = _vPz/_vP;
	      _vRadius = std::sqrt(_vX*_vX+_vY*_vY);
	      _vId     = id;
	    
	      _vNHits++;
	      break;
	    }
	  }
	}
      }
    }

    _nHits = 0;
    _cryEtot = 0.0;

    //some helper variables
    const CaloCluster                       *cluster;
    const std::vector<CaloCrystalHitPtr>    *crystals;
    const CaloCrystalHit                    *crystalHit;
    const RecoCaloDigi                      *recoDigi;
    //    const CaloDigi                    *caloDigi;
    const CaloDigiMC                        *caloDigiMC;
    const SimParticle                       *sim;
    
    for (unsigned int ic=0; ic<caloCrystalHits->size();++ic) {

	   
      CaloCrystalHit const& hit    = caloCrystalHits->at(ic);
      CLHEP::Hep3Vector crystalPos = _calorimeter->crystal(hit.id()).position();

      _cryEtot             += hit.energyDep();
      _cryTime[_nHits]      = hit.time();
      _cryEdep[_nHits]      = hit.energyDep();
	
      _cryPosX[_nHits]      = crystalPos.x()+ 3904;
      _cryPosY[_nHits]      = crystalPos.y();
      _cryPosZ[_nHits]      = crystalPos.z()-10200;
      _cryPosR[_nHits]      = sqrt( _cryPosX[_nHits]*_cryPosX[_nHits] + crystalPos.y()*crystalPos.y() );
      _cryId[_nHits]        = hit.id();
      _crySectionId[_nHits] = _calorimeter->crystal(hit.id()).sectionId();

      int    indexMC, nParticles(0);
      int    isConversion(0);
      
      recoDigi         = hit.recoCaloDigis().at(0).operator ->();

      const CaloDigi	caloDigi  = recoDigi->RODigi();
      indexMC          = caloDigi.index();
      caloDigiMC       = &caloDigiMCCol->at(indexMC);
	
      for (int k=0; k<int(caloDigiMC->nParticles()); ++k){
	sim =   caloDigiMC->simParticle(k).operator ->();
	//	if ( sim->fromGenerator() ){
	const CLHEP::Hep3Vector genPos = sim->startPosition();
	
	if (!_calorimeter->isInsideCalorimeter(genPos)){
	  ++nParticles;
	  GenParticle* gen = (GenParticle*) &(sim->genParticle());
	  if ( gen->generatorId() == GenId::conversionGun ){
	    isConversion = 1;
	  }
	}
      }//end loop on the particles inside the crystalHit

      _cryMCTime    [_nHits] = caloDigiMC->timeFirst();
      _cryMCEdep    [_nHits] = caloDigiMC->totalEDep();
      _cryNParticles[_nHits] = nParticles;
      _cryPsd       [_nHits] = recoDigi->psd();
      _cryIsConv    [_nHits] = isConversion;
      
      ++_nHits;
    }
    
    
    _nCluster = caloClusters->size();

    for (int i=0; i<_nCluster; ++i){
      cluster  = &caloClusters->at(i);

      crystals = &cluster->caloCrystalHitsPtrVector();
      
      int    nCrystals = crystals->size();
      int    indexMC;
      int    isConversion(0);
      
      double   energyMax(0), clusterTime(0), clusterMCMeanTime(0), clusterMCTime(0), eDep, psd;
      
      for (int j=0; j<nCrystals; ++j){
	crystalHit	 = crystals->at(j).operator ->();
	recoDigi         = crystalHit->recoCaloDigis().at(0).operator ->();

	const CaloDigi	caloDigi  = recoDigi->RODigi();
	indexMC          = caloDigi.index();
	caloDigiMC       = &caloDigiMCCol->at(indexMC);
	
	eDep             = crystalHit->energyDep();
	psd              = recoDigi  ->psd();

	if ( (eDep > energyMax) && (psd > _psdThreshold) ){
	  clusterTime       = crystalHit->time();
	  clusterMCMeanTime = caloDigiMC->meanTime();
	  energyMax         = eDep;
	  clusterMCTime     = caloDigiMC->timeFirst();
	}
	
	for (int k=0; k<int(caloDigiMC->nParticles()); ++k){
	  sim =   caloDigiMC->simParticle(k).operator ->();
	  if ( sim->fromGenerator() ){
	    GenParticle* gen = (GenParticle*) &(sim->genParticle());
	    if ( gen->generatorId() == GenId::conversionGun ){
	      isConversion = 1;
	    }
	  }
	}//end loop on the particles inside the crystalHit
      }

      _cluEnergy    [i]     = cluster->energyDep();
      _cluTime      [i]     = clusterTime;
      _cluMCTime    [i]     = clusterMCTime;
      _cluMCMeanTime[i]     = clusterMCMeanTime;
      _cluCrysE     [i]     = energyMax;
      _cluNcrys     [i]     = nCrystals;
      _cluCogX      [i]     = cluster->cog3Vector().x();
      _cluCogY      [i]     = cluster->cog3Vector().y();
      _cluCogR      [i]     = sqrt(_cluCogX[i]*_cluCogX[i] + _cluCogY[i]*_cluCogY[i]);
      _cluCogZ      [i]     = cluster->cog3Vector().z();
      _cluConv      [i]     = isConversion;
      
    }//end filling calo clusters info

    _Ntup->Fill();
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadCaloDigi);


