 //
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
//
// Original author Rob Kutschke
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TDirectory.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2.h"

// Mu2e includes.
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

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



  class CaloCalib : public art::EDAnalyzer {

     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit CaloCalib(fhicl::ParameterSet const& pset);
       virtual ~CaloCalib() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:

       typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
       typedef art::Ptr< CaloHit> CaloHitPtr;
       typedef art::Ptr<SimParticle> SimParticlePtr;


       int _diagLevel;
       int _nProcess;

       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       std::string _caloReadoutModuleLabel;
       std::string _caloCrystalModuleLabel;
       std::string _caloHitMCCrystalPtrLabel;
       std::string _caloClusterModuleLabel;
       std::string _caloClusterAlgorithm;
       std::string _caloClusterSeeding;
       const std::string _producerName;
       std::string _virtualDetectorLabel;
       std::string _trkPatRecModuleLabel;
       std::string _instanceName;

       TrkFitDirection _fdir;


       TH2F* _hviewxy;
       TH2F* _hviewxz;


       TTree* _Ntup;


       int   _evt,_run;

       int   _nGen,_genPdgId[16384],_genCrCode[16384];
       float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];

       int   _nHits,_cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
       float _cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];

       int   _nSim,_motId[16384],_motPdgId[16384],_motcrCode[16384],_motGenIdx[16384];
       float _motmom[16384],_motStartX[16384],_motStartY[16384],_motStartZ[16384],_motStartT[16384];
       float _motTime[16384],_motEdep[16348],_motPosX[16384],_motPosY[16384],_motPosZ[16384];




  };


  CaloCalib::CaloCalib(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel","g4run")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel","generate")),
    _caloReadoutModuleLabel(pset.get<string>("caloReadoutModuleLabel","CaloReadoutHitsMaker")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel","CaloHitsMaker")),
    _caloHitMCCrystalPtrLabel(pset.get<string>("calorimeterHitMCCrystalPtr","CaloHitMCCrystalPtr")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel","CaloClusterMakerNew")),
    _virtualDetectorLabel(pset.get<string>("virtualDetectorName","virtualdetector")),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel","TPRDownstreameMinus")),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _hviewxy(0),_hviewxz(0),
    _Ntup(0)

  {
    _instanceName = _fdir.name();
  }

  void CaloCalib::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("Calo", "Calo");



    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");

    _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
    _Ntup->Branch("genId",        &_genPdgId,     "genId[nGen]/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode[nGen]/I");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomX[nGen]/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX[nGen]/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY[nGen]/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ[nGen]/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT[nGen]/F");

    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryDose",      &_cryDose ,     "cryDose[nCry]/F");
    _Ntup->Branch("crySimIdx",    &_crySimIdx ,   "crySimIdx[nCry]/I");
    _Ntup->Branch("crySimLen",    &_crySimLen ,   "crySimLen[nCry]/I");

    _Ntup->Branch("nSim",         &_nSim ,        "nSim/I");
    _Ntup->Branch("simId",        &_motId ,       "simId[nSim]/I");
    _Ntup->Branch("simPdgId",     &_motPdgId ,    "simPdgId[nSim]/I");
    _Ntup->Branch("simCrCode",    &_motcrCode ,   "simCrCode[nSim]/I");
    _Ntup->Branch("simMom",       &_motmom ,      "simMom[nSim]/F");
    _Ntup->Branch("simStartX",    &_motStartX ,   "simStartX[nSim]/F");
    _Ntup->Branch("simStartY",    &_motStartY ,   "simStartY[nSim]/F");
    _Ntup->Branch("simStartZ",    &_motStartZ ,   "simStartZ[nSim]/F");
    _Ntup->Branch("simStartT",    &_motStartT ,   "simStartT[nSim]/F");
    _Ntup->Branch("simPosX",      &_motPosX ,     "simPosX[nSim]/F");
    _Ntup->Branch("simPosY",      &_motPosY ,     "simPosY[nSim]/F");
    _Ntup->Branch("simPosZ",      &_motPosZ ,     "simPosZ[nSim]/F");
    _Ntup->Branch("simTime",      &_motTime ,     "simTime[nSim]/F");
    _Ntup->Branch("simEdep",      &_motEdep ,     "simEdep[nSim]/F");
    _Ntup->Branch("simGenIdx",    &_motGenIdx ,   "simGenIdx[nSim]/I");

    _hviewxy = tfs->make<TH2F>("viewxy","StepPoint MC hits dist xy",100,-20.,20.,100,-20.,20.);
    _hviewxz = tfs->make<TH2F>("viewxz","StepPoint MC hits dist xz",230,-10.,220.,100,-20.,20.);


  }



  void CaloCalib::endJob(){
  }




  void CaloCalib::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%1000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;


      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      //Calorimeter const & cal = *(GeomHandle<Calorimeter>());


      //Get generated particles
      art::Handle<GenParticleCollection> gensHandle;
      event.getByLabel(_generatorModuleLabel, gensHandle);
      GenParticleCollection const& genParticles(*gensHandle);

      //Get calorimeter readout hits (2 readout / crystal as of today)
      //art::Handle<CaloHitCollection> caloHitsHandle;
      //event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
      //CaloHitCollection const& caloHits(*caloHitsHandle);

      //Get calo crystal hits (average from readouts)
      //art::Handle<CaloHitCollection> CaloHitsHandle;
      //event.getByLabel(_caloCrystalModuleLabel, CaloHitsHandle);
      //CaloHitCollection const& CaloHits(*CaloHitsHandle);

      //const double CrDensity = 4.9e-6;  // in kg/mm3 to be consistent with volume units!
      //const double CrMass    = CrDensity*cal.caloInfo().crystalVolume();




       //--------------------------  Do generated particles --------------------------------


       _evt = event.id().event();
       _run = event.run();

       _nGen = genParticles.size();
       for (unsigned int i=0; i < genParticles.size(); ++i)
       {
           GenParticle const* gen = &genParticles[i];
           _genPdgId[i]   = gen->pdgId();
           _genCrCode[i]  = gen->generatorId().id();
           _genmomX[i]    = gen->momentum().vect().x();
           _genmomY[i]    = gen->momentum().vect().y();
           _genmomZ[i]    = gen->momentum().vect().z();
           _genStartX[i]  = gen->position().x()+ 3904;
           _genStartY[i]  = gen->position().y();
           _genStartZ[i]  = gen->position().z();
           _genStartT[i]  = gen->time();
       }



       //--------------------------  Do calorimeter hits --------------------------------

       _nHits = _nSim = 0;

/*
       for (unsigned int ic=0; ic<CaloHits.size();++ic)
       {
           CaloHit const& hit            = CaloHits.at(ic);
           int diskId                    = cal.crystal(hit.id()).diskID();
           CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
           CaloHit const& caloHit        = *(hit.readouts().at(0));


           CaloHitSimPartMC const& hitSim = caloHitNavigator.sim(caloHit);
           int nPartInside                = hitSim.simParticles().size();



           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryDose[_nHits]      = hit.energyDep() / CrMass / (CLHEP::joule/CLHEP::kg); //dose
           _cryPosX[_nHits]      = crystalPos.x();
           _cryPosY[_nHits]      = crystalPos.y();
           _cryPosZ[_nHits]      = crystalPos.z();
           _cryId[_nHits]        = hit.id();
           _crySectionId[_nHits] = cal.crystal(hit.id()).diskID();
           _crySimIdx[_nHits]    = _nSim;
           _crySimLen[_nHits]    = nPartInside;

           for (int ip=0; ip<nPartInside;++ip)
           {

             art::Ptr<SimParticle> const& mother = hitSim.simParticles().at(ip);

             art::Ptr<SimParticle> grandMother = mother;
             while (grandMother->hasParent()) grandMother = grandMother->parent();
             GenParticle const* generated = grandMother->genParticle() ? grandMother->genParticle().get() : 0;

             CLHEP::Hep3Vector hitSimPos = cal.geomUtil().mu2eToDiskFF(diskId,hitSim.position().at(ip)); //in disk FF frame

             _motId[_nSim]      = mother->id().asInt();
             _motPdgId[_nSim]   = mother->pdgId();
             _motmom[_nSim]     = hitSim.momentum().at(ip);
             _motcrCode[_nSim]  = mother->creationCode();
             _motStartX[_nSim]  = mother->startPosition().x() + 3904; //in Mu2e frame
             _motStartY[_nSim]  = mother->startPosition().y();
             _motStartZ[_nSim]  = mother->startPosition().z() - 10200;  if this is uncommented, the hardwired 10200 should be changed to parametrized version
             _motStartT[_nSim]  = mother->startGlobalTime();
             _motPosX[_nSim]    = hitSimPos.x(); // in disk FF frame
             _motPosY[_nSim]    = hitSimPos.y();
             _motPosZ[_nSim]    = hitSimPos.z();
             _motTime[_nSim]    = hitSim.time().at(ip);
             _motEdep[_nSim]    = hitSim.energyDep().at(ip);

             _motGenIdx[_nSim]  = -1;
             if (generated) _motGenIdx[_nSim] = generated - &(genParticles.at(0));
             ++_nSim;

           }

           ++_nHits;
       }
*/


          _Ntup->Fill();





  }



}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloCalib)


