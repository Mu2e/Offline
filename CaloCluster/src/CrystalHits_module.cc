//
// $Id: CrystalHits_module.cc,v 1.4 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
//
//Original author Giovanni Onorato

// Mu2e includes.
#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

//Root includes
#include "TFile.h"
#include "TTree.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/DataViewImpl.h"


// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "MCDataProducts/inc/StatusG4.hh"
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {


  class CrystalHits : public art::EDAnalyzer {
  public:
    explicit CrystalHits(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),   
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _Ntup(0)
    {
    }
    virtual ~CrystalHits() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:
    int _diagLevel;
    //Label for the generator
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    std::string _clusterSeed;
    std::string _clusterAlgorithm;
    std::string _producerName;

    // Label of the calo readout hits maker
    std::string _caloReadoutModuleLabel;

    // Label of the calo crystal hists maker
    std::string _caloCrystalModuleLabel;

    TTree* _Ntup;

    std::unique_ptr<MCCaloUtilities> CaloManager;


    Int_t //_clNo,
    //        _nCryCl,
    _cryId,
      _vane,
      _crySeedsPdgId[10000],
      _cryPdgIdFirst,
      _cryIsGen,
      _crynStepsFromGen,
      _crySeedsTrkId[10000],
      _cryIsConv[10000],

      _crytrkIdFirst,
      _crynSteps;

    Float_t _evt,
      _cryRow,
      _cryColumn,
      _cryEdep,
      _cryEdepTot,
      _cryT,
      _crySeedsTime[10000],
      _crySeedsEnergy[10000];

    Int_t _SPpdgId[10000],
      _SPIsGen[10000],
      _SPTrkId[10000];

    Float_t _SPx[10000],
      _SPy[10000],
      _SPz[10000],
      _SPu[10000],
      _SPv[10000],
      _SPw[10000],
      _SPT[10000],
      _SPLength[10000],
      _SPE[10000],
      _SPpx[10000],
      _SPpy[10000],
      _SPpz[10000],
      _SPpCosTh[10000],
      _SPpPhi[10000];

  };

  void CrystalHits::beginJob( ) {

  }

  void CrystalHits::analyze(art::Event const& evt ) {

    static int ncalls(0);
    ++ncalls;

    if (ncalls == 1) {

      // cout << "This should be done only in the first event" << endl;


      art::ServiceHandle<art::TFileService> tfs;

      _Ntup        = tfs->make<TTree>("crystalsTree", "crystal hits info");

      _Ntup->Branch("evt", &_evt , "evt/F");

      _Ntup->Branch("cryId", &_cryId, "cryId/I");
      _Ntup->Branch("cryRow", &_cryRow, "cryRow/F");
      _Ntup->Branch("cryColumn", &_cryColumn, "cryColumn/F");
      _Ntup->Branch("vane",&_vane , "vane/I");
      _Ntup->Branch("crynStepsFromGen",&_crynStepsFromGen , "crynStepsFromGen/I");
      _Ntup->Branch("cryEdep",&_cryEdep , "cryEdep/F");
      _Ntup->Branch("cryEdepTot",&_cryEdepTot , "cryEdepTot/F");
      _Ntup->Branch("cryT",&_cryT , "cryT/F");
      _Ntup->Branch("crySeedsPdgId[crynStepsFromGen]",_crySeedsPdgId , "crySeedsPdgId[crynStepsFromGen]/I");
      _Ntup->Branch("crySeedsTrkId[crynStepsFromGen]",_crySeedsTrkId, "crySeedsTrkId[crynStepsFromGen]/I");

      _Ntup->Branch("cryIsConv[crynStepsFromGen]",_cryIsConv, "cryIsConv[crynStepsFromGen]/I");//_crySeedsTime

      _Ntup->Branch("crySeedsTime[crynStepsFromGen]",_crySeedsTime, "crySeedsTime[crynStepsFromGen]/F");//_crySeedsTime
      _Ntup->Branch("crySeedsEnergy[crynStepsFromGen]",_crySeedsEnergy, "crySeedsEnergy[crynStepsFromGen]/F");//_cryIsConv
      _Ntup->Branch("cryIsGen",&_cryIsGen , "cryIsGen/I");
      _Ntup->Branch("cryPdgIdFirst",&_cryPdgIdFirst , "cryPdgIdFirst/I");
      _Ntup->Branch("crytrkIdFirst",&_crytrkIdFirst , "crytrkIdFirst/I");
      _Ntup->Branch("crynSteps",&_crynSteps , "crynSteps/I");

      _Ntup->Branch("SPx[crynSteps]",_SPx , "SPx[crynSteps]/F");
      _Ntup->Branch("SPy[crynSteps]",_SPy , "SPy[crynSteps]/F");
      _Ntup->Branch("SPz[crynSteps]",_SPz , "SPz[crynSteps]/F");
      _Ntup->Branch("SPu[crynSteps]",_SPu , "SPu[crynSteps]/F");
      _Ntup->Branch("SPv[crynSteps]",_SPv , "SPv[crynSteps]/F");
      _Ntup->Branch("SPw[crynSteps]",_SPw , "SPw[crynSteps]/F");
      _Ntup->Branch("SPT[crynSteps]",_SPT , "SPT[crynSteps]/F");
      _Ntup->Branch("SPpdgId[crynSteps]",_SPpdgId , "SPpdgId[crynSteps]/I");
      _Ntup->Branch("SPLength[crynSteps]",_SPLength , "SPLength[crynSteps]/F");
      _Ntup->Branch("SPE[crynSteps]",_SPE , "SPE[crynSteps]/F");
      _Ntup->Branch("SPpx[crynSteps]",_SPpx , "SPpx[crynSteps]/F");
      _Ntup->Branch("SPpy[crynSteps]",_SPpy , "SPpy[crynSteps]/F");
      _Ntup->Branch("SPpz[crynSteps]",_SPpz , "SPpz[crynSteps]/F");
      _Ntup->Branch("SPpCosTh[crynSteps]",_SPpCosTh , "SPpCosTh[crynSteps]/F");
      _Ntup->Branch("SPpPhi[crynSteps]",_SPpPhi , "SPpPhi[crynSteps]/F");
      _Ntup->Branch("SPIsGen[crynSteps]",_SPIsGen , "SPIsGen[crynSteps]/I");
      _Ntup->Branch("SPTrkId[crynSteps]",_SPTrkId , "SPTrkId[crynSteps]/I");

    }

    art::Handle<CaloHitCollection> caloHits;
    evt.getByLabel(_caloReadoutModuleLabel, caloHits);

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
    evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    GeomHandle<VaneCalorimeter> cg;
    if (!( caloHits.isValid())) {
      return;
    }

    if (!caloCrystalHits.isValid()) {
      cout << "NO CaloCrystalHits" << endl;
      return;
    }

    //cout<<"caloCrystalHits->size : "<<caloCrystalHits->size()<<endl;
    for(size_t i=0; i<caloCrystalHits->size(); ++i){
      
      if( _diagLevel>0){
	cout<<"caloCrystalHits->at : "<<i<<endl;
      }
      
      CaloCrystalHit const& cry = (caloCrystalHits->at(i));
      
      if( _diagLevel>0){
	cout<<"preso cry..."<<endl;
      }
      
      _evt = evt.id().event();
      _cryT = cry.time();
      _cryEdep = cry.energyDep();
      _cryEdepTot = cry.energyDepTotal();

      if(cry.readouts().size()==0) {
	if( _diagLevel>0){
	  cout<<"cry.readouts().size()==0"<<endl;
	}
	continue;
      }
      CaloHit const & RO = *(cry.readouts().at(0));
     if( _diagLevel>0){
       cout<<"presa Ro..."<<endl;
     }
      _cryId = cg->crystalByRO(RO.id());
      _vane = cg->vaneByRO(RO.id());
      
      PtrStepPointMCVector const & mcptr(mcptrHandle->at(cry.readouts().at(0).key()));
      if( _diagLevel>0){
	cout<<"preso mcptr..."<<endl;
      }
      
      if(mcptr.size()==0) continue;

      _crynSteps = mcptr.size();

      if( _diagLevel>0){
	cout<<"-------------- 1 -------------"<<endl
	    << "mcptr.size() = "<< mcptr.size()<<endl
	    <<"-------------- 2 -------------"<<endl;
      }
      float earliest = 100000;

      if( _diagLevel>0){
	cout<<"_crynSteps = "<<_crynSteps<<endl;
      }

      _crynStepsFromGen = 0;
      std::map<unsigned int, unsigned int> searchSeeds;//the key is the trkId and the value is the PdgId

      for ( int k = 0; k < _crynSteps; ++k) {
	const StepPointMC & SP = *mcptr[k];
	if(_diagLevel>0){
	  cout<<"step point at "<<k<<endl;
	  SP.print();
	}


	art::Ptr<SimParticle> const& simp = SP.simParticle();
	SimParticle const& sim = *simp;
	
	 if( _diagLevel>0){
	   cout<<"simParticle pointer taken..."<<endl
	       <<" pdgId = "<< sim.pdgId() <<endl;
	 }
	 
	if(sim.fromGenerator()){
	  searchSeeds[SP.trackId().asUint()] = k;//simp->pdgId();
	}
      }

      _crynStepsFromGen = searchSeeds.size();
      if( _diagLevel>0){
	cout<<"crynStepsFromGen = "<<_crynStepsFromGen<<endl;
      }

      size_t count = 0;
      for(std::map<unsigned int, unsigned int> ::iterator it = searchSeeds.begin(); it != searchSeeds.end(); ++it){
	unsigned int index = it->second;
	const StepPointMC & SP = *(mcptr.at(index));
	art::Ptr<SimParticle> const& simp = SP.simParticle();
	SimParticle const& sim = *simp;

	_cryIsConv[count] = 999999;
	
	if(sim.genParticle().isNonnull()){
	  GenParticle const& gen =*sim.genParticle();
	  GenId genId = gen.generatorId();
	  
	  if(genId==GenId::conversionGun){
	    _cryIsConv[count] = 1;
	  }else {
	    _cryIsConv[count] = 0;
	  }
	}

	_crySeedsEnergy[count] = SP.momentum().mag();
	_crySeedsTime[count] = SP.time();
	_crySeedsPdgId[count] = sim.pdgId();
	_crySeedsTrkId[count] = it->first;
	++count;
      }


      for ( int k = 0; k < _crynSteps; ++k) {
	if( _diagLevel>0){
	  cout<<"crynSteps->at : "<<k<<endl;
	}
	const StepPointMC & SP = *(mcptr.at(k));
	art::Ptr<SimParticle> const& simptr = SP.simParticle();
	SimParticle const& sim = *simptr;

	if( _diagLevel>0){
	  cout<<"srtep point preso..."<<endl;
	}
	_SPx[k] = SP.position().x();
	_SPy[k] = SP.position().y();
	_SPz[k] = SP.position().z();
	CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame( _cryId, SP.position());
	_SPu[k] = cryFrame.x();
	_SPv[k] = cryFrame.y();
	_SPw[k] = cryFrame.z();
	_SPLength[k] = SP.stepLength();
	_SPE[k] = SP.totalEDep();
	_SPpx[k] = SP.momentum().x();
	_SPpy[k] = SP.momentum().y();
	_SPpz[k] = SP.momentum().z();
	_SPpCosTh[k] = SP.momentum().cosTheta();
	_SPpPhi[k] = SP.momentum().phi();
	_SPT[k] = SP.time();
	_SPpdgId[k] = sim.pdgId();
	_SPIsGen[k] = sim.fromGenerator();
	_SPTrkId[k] = SP.trackId().asInt();
	if (_SPT[k] < earliest) {
	  earliest  = _SPT[k];
	  _cryPdgIdFirst = _SPpdgId[k];
	  _cryIsGen = _SPIsGen[k];
	  _crytrkIdFirst = _SPTrkId[k];
	}

      }

      _Ntup->Fill();

    }
    cout << "Event "<<evt.id().event()<<" CrystalHits done..."<<endl;
  } // end of analyze

  void CrystalHits::endJob() {
  }

}

using mu2e::CrystalHits;
DEFINE_ART_MODULE(CrystalHits);

