//
// A module to evaluate the normalization of background to simulate
//
// $Id: BkgNorm_module.cc,v 1.8 2012/02/06 23:56:32 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/02/06 23:56:32 $
//
// Original author Gianni Onorato
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {


  class BkgNorm : public art::EDAnalyzer {
  public:
    explicit BkgNorm(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _tNtup(0),
      _cNtup(0),
      _nDevices(36),
      _nSectors(6),
      _nLayers(2),
      _nStrawsPerLay(50),
      _nBadG4Status(0),
      _nOverflow(0),
      _nKilled(0),
      _totalcputime(0),
      _totalrealtime(0)
    {
    }
    virtual ~BkgNorm() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    void doTracker(art::Event const& evt, bool skip);
    void doCalorimeter(art::Event const& evt, bool skip);
    void doStoppingTarget(art::Event const& evt);

    // Diagnostic level
    int _diagLevel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    // Label of the Calorimeter modules;
    std::string _caloReadoutModuleLabel;
    std::string _caloCrystalModuleLabel;


    double _minimumEnergy; //minimum energy deposition of hits

    TNtuple* _tNtup, *_cNtup;
    const int _nDevices, _nSectors, _nLayers, _nStrawsPerLay;

    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;

    bool _skipEvent;

  };


  void BkgNorm::beginJob( ) {
  }

  void BkgNorm::analyze(art::Event const& evt ) {

    static int ncalls(0);
    ++ncalls;

    art::Handle<StatusG4> g4StatusHandle;
    evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting BkgNorm::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    if (g4Status.overflowSimParticles()) {
      ++_nOverflow;
      mf::LogError("G4")
        << "Aborting BkgNorm::analyze due to overflow of particles\n"
        << g4Status;
      return;
    }

    if (g4Status.nKilledStepLimit() > 0) {
      ++_nKilled;
      mf::LogError("G4")
        << "Aborting BkgNorm::analyze due to nkilledStepLimit reached\n"
        << g4Status;
      return;
    }

    _totalcputime += g4Status.cpuTime();
    _totalrealtime += g4Status.realTime();

    if (ncalls == 1) {

      //      art::Handle<PhysicalVolumeInfoCollection> volumes;
      // evt.getRun().getByLabel(_g4ModuleLabel, volumes);

      // for (size_t i=0; i < volumes->size(); ++i) {

      //        PhysicalVolumeInfo const& volInfo = volumes->at(i);
      //  cout << i << '\t' <<volInfo.name() << volInfo.copyNo() << endl;
      // }



      art::ServiceHandle<art::TFileService> tfs;

      // "evt:run:"; //Event info (2 entries)
      // "time:dt:eDep:lay:dev:sec:strawId:strawX:strawY:strawZ:"; //StrawHit info (10 entries)
      // "trkPdgId:trkP:trkIsGen:trkStartVolume:trkStepFromEva:EvaIsGen"; //Track making hits info (6 entries)
      // "genPdgId:genId:genP:genE:genX:genY:genZ:genT:genPhi:genCosth:"; //Generated particle info (10 entries)
      // "dau1PdgId:dau1P:dau1StartVolume:"; //First Daughter of generated particle info (3 entries)

      _tNtup        = tfs->make<TNtuple>( "StrawHits", "Straw Ntuple", "evt:run:time:dt:eDep:lay:dev:sec:strawId:strawX:strawY:strawZ:trkPdgId:trkP:trkIsGen:trkStartVolume:trkStepFromEva:EvaIsGen:genPdgId:genId:genP:genE:genX:genY:genZ:genT:genPhi:genCosth:dau1PdgId:dau1P:dau1StartVolume");
      _cNtup        = tfs->make<TNtuple>( "CaloHits", "calo Ntupla", "evt:run:time:eDep:vane:crId:trkPdgId:trkP:trkIsGen:trkStartVolume:trkStepFromEva:EvaIsGen:genPdgId:genId:genP:genE:genX:genY:genZ:genT:genPhi:genCosth:dau1PdgId:dau1P:dau1StartVolume");
   }

    doStoppingTarget(evt);

    doTracker(evt, _skipEvent);
    doCalorimeter(evt, _skipEvent);

  } // end of analyze

  void BkgNorm::endJob() {
    cout << "BkgNorm::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
         << "\nBkgNorm::endJob Number of overflow events "
         << "due to too many particles in G4: "
         << _nOverflow
         << "\nBkgNorm::endJob Number of events with killed particles "
         << "due to too many steps in G4: "
         << _nKilled
         << "\nBkgNorm::endJob total CpuTime "
         << _totalcputime
         << "\nBkgNorm::endJob total RealTime "
         << _totalrealtime
         << endl;
  }


  void BkgNorm::doTracker(art::Event const& evt, bool skip) {

    if (skip) return;

    const Tracker& tracker = getTrackerOrThrow();

    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    if (!(hits->size() == hits_truth->size() &&
          hits_mcptr->size() == hits->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits: " << hits->size()
        << " MCTruthStrawHits: " << hits_truth->size()
        << " MCPtr: " << hits_mcptr->size();
    }

    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);

    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    size_t nStrawPerEvent = hits->size();

    for (size_t i=0; i<nStrawPerEvent; ++i) {

      // Access data
      StrawHit             const&      hit(hits->at(i));
      PtrStepPointMCVector const&    mcptr(hits_mcptr->at(i));

      double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      if (hitEnergy < _minimumEnergy) continue;

      //Get hit straw
      StrawIndex si = hit.strawIndex();
      Straw str = tracker.getStraw(si);
      StrawId sid = str.id();
      LayerId lid = sid.getLayerId();
      DeviceId did = sid.getDeviceId();
      SectorId secid = sid.getSectorId();

      //Get coordinates of the hit:

      //X, Y and Z coordinate of the straw middle point

      const CLHEP::Hep3Vector stMidPoint3 = str.getMidPoint();
      double xc = stMidPoint3.getX();
      double yc = stMidPoint3.getY();
      double z = stMidPoint3.getZ();

      //time of the hit
      double hitTime = hit.time();

      float tntpArray[31];
      int idx(0);
      tntpArray[idx++] = evt.id().event();
      tntpArray[idx++] = evt.run();
      tntpArray[idx++] = hitTime;
      tntpArray[idx++] = hit.dt();
      tntpArray[idx++] = hitEnergy;
      tntpArray[idx++] = lid.getLayer();
      tntpArray[idx++] = did;
      tntpArray[idx++] = secid.getSector();
      tntpArray[idx++] = sid.getStraw();
      tntpArray[idx++] = xc;
      tntpArray[idx++] = yc;
      tntpArray[idx++] = z;

      bool notFirstTrack = true;

      size_t j = 0;

      CLHEP::Hep3Vector const& strDir = str.direction();

      double strRadius = str.getRadius();

      while (notFirstTrack && j<mcptr.size()) {
        if (j==mcptr.size()-1) {
          j=0;
          break;
        }
        StepPointMC const& mchit = *mcptr[j];
        SimParticle const& sim = simParticles->at(mchit.trackId());
        CLHEP::Hep3Vector const& StartPos = sim.startPosition();
        LinePointPCA lppca(stMidPoint3, strDir, StartPos);
        double insideDistance = lppca.dca();
        if (insideDistance >= strRadius) {
          notFirstTrack = false;
          break;
        } else {
          ++j;
        }
      }

      StepPointMC const& mchit = *mcptr[j];
      SimParticle const& sim = simParticles->at(mchit.trackId());
      tntpArray[idx++] = sim.pdgId();
      tntpArray[idx++] = sim.startMomentum().vect().mag();
      tntpArray[idx++] = sim.fromGenerator();
      tntpArray[idx++] = sim.startVolumeIndex();
      int nEvolutionSteps = 0;
      SimParticleCollection::key_type Dau1Idx = SimParticleCollection::key_type(0);
      if (!sim.fromGenerator()) {
        bool notEva = true;
        SimParticle& baby = const_cast<SimParticle&>(sim);
        while (notEva) {
          if (!baby.hasParent()) {
            tntpArray[idx++] = nEvolutionSteps;
            tntpArray[idx++] = 0;
            notEva = false;
            break;
          }
          SimParticle & mommy = const_cast<SimParticle&>(*baby.parent());
          nEvolutionSteps++;
          if (mommy.fromGenerator()) {
            tntpArray[idx++] = nEvolutionSteps;
            tntpArray[idx++] = 1;
            notEva = false;
            Dau1Idx = baby.id();
            break;
          } else {
            baby = mommy;
          }
        }
      } else {
        tntpArray[idx++] = nEvolutionSteps;
        tntpArray[idx++] = 1;
      }


      SimParticleCollection::key_type idxInSim = SimParticleCollection::key_type(1);
      SimParticle const& geninSim = simParticles->at(idxInSim);
      GenParticle const& gen = genParticles->at(geninSim.generatorIndex());
      tntpArray[idx++] = gen.pdgId();
      tntpArray[idx++] = gen.generatorId().id();
      tntpArray[idx++] = gen.momentum().vect().mag();
      tntpArray[idx++] = gen.momentum().e();
      tntpArray[idx++] = gen.position().x();
      tntpArray[idx++] = gen.position().y();
      tntpArray[idx++] = gen.position().z();
      tntpArray[idx++] = gen.time();
      tntpArray[idx++] = gen.momentum().cosTheta();
      tntpArray[idx++] = gen.momentum().phi();

      if (Dau1Idx != SimParticleCollection::key_type(0)) {
        SimParticle const& Dau1 = simParticles->at(Dau1Idx);
        tntpArray[idx++] = Dau1.pdgId();
        tntpArray[idx++] = Dau1.startMomentum().vect().mag();
        tntpArray[idx++] = Dau1.startVolumeIndex();
      } else {
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
      }

      _tNtup->Fill(tntpArray);


    } //end of Strawhits loop

  } // end of doTracker

  void BkgNorm::doCalorimeter(art::Event const& evt, bool skip) {

    if (skip) return;

    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    art::Handle<CaloHitCollection> caloHits;
    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;

    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);
    evt.getByLabel(_caloReadoutModuleLabel, caloHits);
    evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
    if (!( caloHits.isValid())) {
      return;
    }

    if (!caloCrystalHits.isValid()) {
      cout << "NO CaloCrystalHits" << endl;
      return;
    }

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );
    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }


    if (caloCrystalHits->size()>0) {

      for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {

        CaloCrystalHit const & hit = (*caloCrystalHits).at(i);

        std::vector<art::Ptr<CaloHit> > const & ROIds  = hit.readouts();

        if (hit.energyDep() < _minimumEnergy) continue;
        if (ROIds.size() < 1) continue;

        bool readCryOnce(false);
        float cntpArray[25];
        int idx(0);

        double firstHitTime = 100000;
        CLHEP::Hep3Vector firstHitPos(0,0,0);
        CLHEP::Hep3Vector cryFrame(0,0,0);
        size_t firstTrackIndex = 0;
        size_t fTCollPos = 0;

        for (size_t it = 0;
             it < ROIds.size() ; ++it ) {

          size_t collectionPosition = ROIds.at(it).key();
          CaloHit const & thehit = *ROIds.at(it);

          if (!readCryOnce) {
            CLHEP::Hep3Vector cryCenter =  cg->getCrystalOriginByRO(thehit.id());
            int vane = cg->getVaneByRO(thehit.id());
            cntpArray[idx++] = evt.id().event();
            cntpArray[idx++] = evt.run();
            cntpArray[idx++] = hit.time();
            cntpArray[idx++] = hit.energyDep();
            cntpArray[idx++] = vane;
            cntpArray[idx++] = cg->getCrystalByRO(thehit.id());
            readCryOnce = true;
          }

          PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
          size_t nHitsPerCrystal = mcptr.size();

          for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

            StepPointMC const& mchit = *mcptr[j2];
            if (mchit.time() < firstHitTime) {
              firstHitTime = mchit.time();
              firstTrackIndex = j2;
              fTCollPos = collectionPosition;
            }
          }
        }


        // The simulated particle that made this hit.
        PtrStepPointMCVector const & mcptr(hits_mcptr->at(fTCollPos));
        StepPointMC const& mchit = *mcptr[firstTrackIndex];
        SimParticleCollection::key_type trackId(mchit.trackId());
        SimParticle const& sim = simParticles->at(trackId);
        cntpArray[idx++] = sim.pdgId();
        cntpArray[idx++] = sim.startMomentum().vect().mag();
        cntpArray[idx++] = sim.fromGenerator();
        cntpArray[idx++] = sim.startVolumeIndex();
        int nEvolutionSteps = 0;
        SimParticleCollection::key_type Dau1Idx = SimParticleCollection::key_type(0);
        if (!sim.fromGenerator()) {
        bool notEva = true;
        SimParticle& baby = const_cast<SimParticle&>(sim);
        while (notEva) {
          if (!baby.hasParent()) {
            cntpArray[idx++] = nEvolutionSteps;
            cntpArray[idx++] = 0;
            notEva = false;
            break;
          }
          SimParticle & mommy = const_cast<SimParticle&>(*baby.parent());
          nEvolutionSteps++;
          if (mommy.fromGenerator()) {
            cntpArray[idx++] = nEvolutionSteps;
            cntpArray[idx++] = 1;
            notEva = false;
            Dau1Idx = baby.id();
            break;
          } else {
            baby = mommy;
          }
        }
        } else {
          cntpArray[idx++] = nEvolutionSteps;
          cntpArray[idx++] = 1;
        }

        SimParticleCollection::key_type idxInSim = SimParticleCollection::key_type(1);
        SimParticle const& geninSim = simParticles->at(idxInSim);
        GenParticle const& gen = genParticles->at(geninSim.generatorIndex());
        cntpArray[idx++] = gen.pdgId();
        cntpArray[idx++] = gen.generatorId().id();
        cntpArray[idx++] = gen.momentum().vect().mag();
        cntpArray[idx++] = gen.momentum().e();
        cntpArray[idx++] = gen.position().x();
        cntpArray[idx++] = gen.position().y();
        cntpArray[idx++] = gen.position().z();
        cntpArray[idx++] = gen.time();
        cntpArray[idx++] = gen.momentum().cosTheta();
        cntpArray[idx++] = gen.momentum().phi();

        if (Dau1Idx != SimParticleCollection::key_type(0)) {
          SimParticle const& Dau1 = simParticles->at(Dau1Idx);
          cntpArray[idx++] = Dau1.pdgId();
          cntpArray[idx++] = Dau1.startMomentum().vect().mag();
          cntpArray[idx++] = Dau1.startVolumeIndex();
        } else {
          cntpArray[idx++] = 0;
          cntpArray[idx++] = 0;
          cntpArray[idx++] = 0;
        }

        _cNtup->Fill(cntpArray);

      }
    }

  }


  void BkgNorm::doStoppingTarget(const art::Event& event) {

    bool generatedStopped = false;

    // Find original G4 steps in the stopping target
    art::Handle<StepPointMCCollection> sthits;
    event.getByLabel(_g4ModuleLabel,"stoppingtarget",sthits);

    // SimParticles container
    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    if( !(simParticles.isValid()) || simParticles->empty() ) return;

    art::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByLabel(_g4ModuleLabel, volumes);

    set<SimParticleCollection::key_type> stoppedtracks;

    // Loop over all hits in the stopping target
    for ( size_t i=0; i<sthits->size(); ++i ){

      // This is G4 hit (step) in the target
      const StepPointMC& hit = (*sthits)[i];

      SimParticleCollection::key_type trackId = hit.trackId();

      SimParticle const* sim = simParticles->getOrNull(trackId);
      if( !sim ) continue;

      PhysicalVolumeInfo const& volInfo = volumes->at(sim->endVolumeIndex());

      if ( sim->fromGenerator() && (sim->pdgId() == 13 || sim->pdgId() == -13)) {
        if ( volInfo.name() == "TargetFoil_" ) {
        generatedStopped = true;
        }
      }
    }
    _skipEvent = generatedStopped;
  }  // end doStoppingTarget

}

using mu2e::BkgNorm;
DEFINE_ART_MODULE(BkgNorm)
