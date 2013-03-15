//
// $Id: ClusterParticleTrajectory_module.cc,v 1.5 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
//
//Original author Giovanni Onoratto

// Mu2e includes.
#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
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


  class ClusterParticleTrajectory : public art::EDAnalyzer {
  public:
    explicit ClusterParticleTrajectory(fhicl::ParameterSet const& pset):
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _clusterModuleLabel(pset.get<std::string>("clusterModuleLabel", "makeCaloCluster")),
      _clusterSeed(pset.get<std::string>("caloClusterSeeding", "energy")),
      _clusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
      _producerName("Algo"+TOUpper(_clusterAlgorithm)+"SeededBy"+TOUpper(_clusterSeed)),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _Ntup(0)
    {
    }
    virtual ~ClusterParticleTrajectory() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    std::string _clusterModuleLabel;
    std::string _clusterSeed;
    std::string _clusterAlgorithm;
    std::string _producerName;

    // Label of the calo readout hits maker
    std::string _caloReadoutModuleLabel;

    // Label of the calo crystal hists maker
    std::string _caloCrystalModuleLabel;

    TTree* _Ntup;

    std::unique_ptr<MCCaloUtilities> CaloManager;


    Int_t _clNo,
      _nCryCl,
      _cryId,
      _vane,
      _cryPdgId,
      _cryIsGen,
      _crytrkId,
      _crynSteps;

    Float_t _evt,
        _clE,
        _clT,
        _clCOGx,
        _clCOGy,
        _clCOGz,
      _clCryEmaxRow,
      _clCryEmaxColumn,
      _clSize,
      _cryRow,
      _cryColumn,
        _cryEdep,
        _cryEdepTot,
        _cryT;

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

  void ClusterParticleTrajectory::beginJob( ) {

  }

  void ClusterParticleTrajectory::analyze(art::Event const& evt ) {

    static int ncalls(0);
    ++ncalls;

    if (ncalls == 1) {

      // cout << "This should be done only in the first event" << endl;


      art::ServiceHandle<art::TFileService> tfs;

      _Ntup        = tfs->make<TTree>("ClusterTrj", "Cluster trajectory info");

      _Ntup->Branch("evt", &_evt , "evt/F");
      _Ntup->Branch("clNo",&_clNo , "clNo/I");
      _Ntup->Branch("clE",&_clE , "clE/F");
      _Ntup->Branch("clT",&_clT , "clT/F");
      _Ntup->Branch("clCOGx",&_clCOGx , "clCOGx/F");
      _Ntup->Branch("clCOGy",&_clCOGy , "clCOGy/F");
      _Ntup->Branch("clCOGz",&_clCOGz , "clCOGz/F");
      _Ntup->Branch("clCryEmaxRow",&_clCryEmaxRow , "clCryEmaxRow/F");
      _Ntup->Branch("clCryEmaxColumn",&_clCryEmaxColumn , "clCryEmaxColumn/F");
      _Ntup->Branch("cryRow",&_cryRow , "cryRow/F");
      _Ntup->Branch("cryColumn",&_cryColumn , "cryColumn/F");
      _Ntup->Branch("clSize",&_clSize , "clSize/F");

      _Ntup->Branch("nCryCl",&_nCryCl , "nCryCl/I");
      _Ntup->Branch("cryId", &_cryId, "cryId/I");
      _Ntup->Branch("vane",&_vane , "vane/I");
      _Ntup->Branch("cryEdep",&_cryEdep , "cryEdep/F");
      _Ntup->Branch("cryEdepTot",&_cryEdepTot , "cryEdepTot/F");
      _Ntup->Branch("cryT",&_cryT , "cryT/F");
      _Ntup->Branch("cryPdgId",&_cryPdgId , "cryPdgId/I");
      _Ntup->Branch("cryIsGen",&_cryIsGen , "cryIsGen/I");
      _Ntup->Branch("crytrkId",&_crytrkId , "crytrkId/I");
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


    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
    evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

    art::Handle<CaloClusterCollection> caloClusters;
    evt.getByLabel(_clusterModuleLabel,_producerName,caloClusters );

    GeomHandle<VaneCalorimeter> cg;

    _evt = evt.id().event();

    for (size_t i = 0; i < caloClusters->size() ; ++i ) {

      CaloCluster const& cl = caloClusters->at(i);
      CaloClusterTools cluTool(cl);

      _clNo = i;
      _clE = cl.energyDep();
      _clT = cl.time();
      _clCOGx = cl.cog3Vector().x();
      _clCOGy = cl.cog3Vector().y();
      _clCOGz = cl.cog3Vector().z();
      _clCryEmaxRow = cluTool.cryEnergydepMaxRow() ;
      _clCryEmaxColumn = cluTool.cryEnergydepMaxColumn() ;
      _clSize = cl.size();
      _nCryCl = cl.size();

      CaloCrystalHitPtrVector CryPtrVec = cl.caloCrystalHitsPtrVector();

      for (int j = 0; j < _nCryCl; ++j) {

        CaloCrystalHit const& cry = *(CryPtrVec.at(j));

        _cryT = cry.time();
        _cryEdep = cry.energyDep();
        _cryEdepTot = cry.energyDepTotal();

        CaloHit const & RO = *(cry.readouts().at(0));

        _cryId = cg->crystalByRO(RO.id());
        _vane = cg->vaneByRO(RO.id());
	_cryRow = cg->crystalRByRO(RO.id()) ;
	_cryColumn = cg->crystalZByRO(RO.id());

        PtrStepPointMCVector const & mcptr(mcptrHandle->at(cry.readouts().at(0).key()));

        _crynSteps = mcptr.size();

//        cout<<"-------------- 1 -------------"<<endl;
//        cout<< "mcptr.size() = "<< mcptr.size()<<endl;
//        cout<<"-------------- 2 -------------"<<endl;

        float earliest = 100000;

        for ( int k = 0; k < _crynSteps; ++k) {

          StepPointMC const& SP = *(mcptr.at(k));

          _SPx[k] = SP.position().x();
          _SPy[k] = SP.position().y();
          _SPz[k] = SP.position().z();
          CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(RO.id(), SP.position());
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
          SimParticle const& sim = *(simParticles->getOrNull(SP.trackId()));
          _SPpdgId[k] = sim.pdgId();
          _SPIsGen[k] = sim.fromGenerator();
          _SPTrkId[k] = SP.trackId().asInt();
          if (_SPT[k] < earliest) {
            _cryPdgId = _SPpdgId[k];
            _cryIsGen = _SPIsGen[k];
            _crytrkId = _SPTrkId[k];
          }

        }

        _Ntup->Fill();

      }
    }
  } // end of analyze

  void ClusterParticleTrajectory::endJob() {
  }

}

using mu2e::ClusterParticleTrajectory;
DEFINE_ART_MODULE(ClusterParticleTrajectory);

