//
// A module to study background rates in the detector subsystems.
//
// $Id: BkgRates_module.cc,v 1.18 2011/07/25 20:51:24 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/07/25 20:51:24 $
//
// Original author Gianni Onorato
//

#include "Analyses/inc/MCCaloUtilities.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
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


  class BkgRates : public art::EDAnalyzer {
  public:
    explicit BkgRates(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _skipStoppedParticle(pset.get<bool>("skipStoppedParticle",false)),
      _nAnalyzed(0),
      _hHitMult(0),
      _hStrawEvt(0),
      _hStrawEvtZ1(0),
      _hStrawEvtZ2(0),
      _hRateUZ(0),
      _hRateU(0),
      _hRateUT(0),
      _hRateMaxU(0),
      _hCaloHitMult(0),
      _hCryEvt(0),
      _hCryEvtZ1(0),
      _hCryEvtZ2(0),
      _tNtup(0),
      _cNtup(0),
      _tgtNtup(0),
      _nDevices(36),
      _nSectors(6),
      _nLayers(2),
      _nStrawsPerLay(50),
      _nVanes(4),
      _nZCryPerVane(44),
      _nRCryPerVane(12),
      _nBadG4Status(0),
      _nOverflow(0),
      _nKilled(0),
      _totalcputime(0),
      _totalrealtime(0)
    {
    }
    virtual ~BkgRates() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    void doTracker(art::Event const& evt, bool skip);
    void doITracker(art::Event const& evt, bool skip);

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

    // Label of the calo readout hits maker
    std::string _caloReadoutModuleLabel;

    // Label of the calo crystal hists maker
    std::string _caloCrystalModuleLabel;

    double _minimumEnergy; //minimum energy deposition of hits
    bool _skipStoppedParticle;

    //number of analyzed events
    int _nAnalyzed;


    TH1F* _hHitMult;
    TH1F* _hStrawEvt;
    TH1F* _hStrawEvtZ1;
    TH1F* _hStrawEvtZ2;
    TH2F* _hRateUZ;
    TH1F* _hRateU;
    TH2F* _hRateUT;
    TH1F* _hRateMaxU;
    TH1F* _hCaloHitMult;
    TH1F* _hCryEvt;
    TH1F* _hCryEvtZ1;
    TH1F* _hCryEvtZ2;

    TNtuple* _tNtup;
    TNtuple* _cNtup;
    TNtuple* _tgtNtup;

    const int _nDevices, _nSectors, _nLayers, _nStrawsPerLay;
    const int _nVanes, _nZCryPerVane, _nRCryPerVane;

    vector<TH2F*> _hStrawRates;
    vector<TH2F*> _hCrystalRates;

    std::auto_ptr<MCCaloUtilities> CaloManager;

    bool _skipEvent;

    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;

  };

  typedef  list< pair <SimParticleCollection::key_type , double> >  PairList;



  bool SortByEnergy(const pair<SimParticleCollection::key_type,double> & a,
                    const pair<SimParticleCollection::key_type,double> & b ) {
    return a.second < b.second;
  }


  void PairListAdd(PairList& list, SimParticleCollection::key_type prt, double e) {

    for (PairList::iterator it = list.begin();
         it != list.end(); ++it) {
      if (it->first == prt) {
        it->second += e;
        return;
      }
    }
    list.push_back(pair<SimParticleCollection::key_type,double>(prt,e));
    return;
  }


  void BkgRates::beginJob( ) {

    CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

    art::ServiceHandle<art::TFileService> tfs;

    _hHitMult       = tfs->make<TH1F>( "hHitMult",    "Multiplicity of g4 hit per Straw ",       100,    0.,  100. );
    _hCaloHitMult   = tfs->make<TH1F>( "hCaloHitMult","Multiplicity of g4 hit per Crystal ",     100,    0.,  100. );
    _hStrawEvt      = tfs->make<TH1F>( "hStrawEvt",   "Multiplicity of straws per event ",       200,    0., 2000. );
    _hStrawEvtZ1    = tfs->make<TH1F>( "hStrawEvtZ1", "Multiplicity of straws per event zoom1",  100,    0.,  500. );
    _hStrawEvtZ2    = tfs->make<TH1F>( "hStrawEvtZ2", "Multiplicity of straws per event zoom2",  100,    0.,  100. );
    _hRateUZ        = tfs->make<TH2F>( "hRateUZ",     "Straw hit in u and z coordinates",         80,  380.,  690., 40, -1550., 1550.);
    _hRateU         = tfs->make<TH1F>( "hRateU",      "Straw hit in u coordinate",               100,  380.,  690. );
    _hRateUT        = tfs->make<TH2F>( "hRateUT",     "Straw hit in u coordinate and time",       80,  380.,  690., 40,   700., 1900.);
    _hCryEvt        = tfs->make<TH1F>( "hCryEvt",     "Multiplicity of Crystal per event",       200,    0., 2000.);
    _hCryEvtZ1      = tfs->make<TH1F>( "hCryEvtZ1",   "Multiplicity of Crystal per event zoom1", 100,    0.,  500.);
    _hCryEvtZ2      = tfs->make<TH1F>( "hCryEvtZ2",   "Multiplicity of Crystal per event zoom2",  50,    0.,   50.);

    for (int i=0; i<_nDevices; ++i) {
      stringstream name, descr;
      name << "hRateDev" << i;
      descr << "Straw rates in device " << i;

      _hStrawRates.push_back(tfs->make<TH2F>(name.str().c_str(),
                                             descr.str().c_str(),
                                             _nStrawsPerLay, 0, _nStrawsPerLay,
                                             _nLayers*_nSectors, 0, _nLayers*_nSectors));
      _hStrawRates[i]->Sumw2();
    }

    for (int i=0; i<_nVanes; ++i) {
      stringstream name, descr;
      name << "hRateVane" << i;
      descr << "Crystal rates in vane " << i;
      _hCrystalRates.push_back(tfs->make<TH2F>(name.str().c_str(),
                                               descr.str().c_str(),
                                               _nZCryPerVane, 0, _nZCryPerVane,
                                               _nRCryPerVane, 0, _nRCryPerVane));
      _hCrystalRates[i]->Sumw2();
    }

  }

  void BkgRates::analyze(art::Event const& evt ) {

    ++_nAnalyzed;

    //*****test code******
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

    art::ServiceHandle<GeometryService> geom;

    if (ncalls == 1) {

      // cout << "This should be done only in the first event" << endl;


      art::ServiceHandle<art::TFileService> tfs;

      if (geom->hasElement<TTracker>()) {
        _tNtup        = tfs->make<TNtuple>( "StrawHits", "Straw Ntuple",
                                            "evt:run:time:dt:eDep:lay:dev:sec:strawId:MChitX:MChitY:v:vMC:z:nTrk:t1trkId:t1pdgId:t1eDep:t1isGen:t1StapFromEva:t1P:t1StartVolume:t2trkId:t2pdgId:t2eDep:t2isGen:t2StapFromEva:t2P:t2StartVolume:t3trkId:t3pdgId:t3eDep:t3isGen:t3StapFromEva:t3P:t3StartVolume:trkId:trkPdgId:trkP:trkIsGen:trkStartVolume:trkStepFromEva:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime:dau1PdgId:dau1P:dauStartVolume:driftTime:driftDist" );
      } else if(geom->hasElement<ITracker>()) {
        _tNtup        = tfs->make<TNtuple>( "CellHits", "Cell Ntuple",
                                            "evt:run:time:eDep:lay:superlay:cellId:MChitX:MChitY:wireZMC:nTrk:t1trkId:t1pdgId:t1en:t1isGen:t2trkId:t2pdgId:t2en:t2isGen:t3trkId:t3pdgId:t3en:t3isGen:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime:driftTime:driftDist" );
      }
      _cNtup        = tfs->make<TNtuple>( "CaloHits", "Calo Ntuple",
                                          "evt:run:crTime:crE:crRad:crId:crVane:crX:crY:crZ:ESwr:EOutVane:NtrkOutside:OutsideE1:OutsidePdg1:OutsideE2:OutsidePdg2:OutsideE3:OutsidePdg3:EOutsideAll:EGen:GenHit1x:GenHit1y:GenHit1z:cryFramex:cryFramey:cryFramez:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime" );
      _tgtNtup      = tfs->make<TNtuple>( "ST", "Particle dead in ST ntuple",
                                          "evt:run:time:x:y:z:isGen:pdgId:trkId:foil");

    }

    doStoppingTarget(evt);

    if (geom->hasElement<ITracker>()) {
      //      cout << "ITracker selected" << endl;
      doITracker(evt, _skipEvent);
    }
    if (geom->hasElement<TTracker>()) {
      //      cout << "TTracker selected" << endl;
      doTracker(evt, _skipEvent);
    }
    doCalorimeter(evt, _skipEvent);



  } // end of analyze

  void BkgRates::endJob() {
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
    
  void BkgRates::doTracker(art::Event const& evt, bool skip) {

    if (skip) return;

    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.

    art::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    /*
    //Get and printout all straw (just to check)
    const std::deque<Straw>& allstraws = tracker.getAllStraws();
    for (size_t i = 0; i< allstraws.size(); ++i) {
    Straw str = allstraws[i];
    cout << "Index " << str.Index()
    << "\t" << str.id() << endl;
    }
    */

    if (!(hits->size() == hits_truth->size() &&
          hits_mcptr->size() == hits->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits: " << hits->size()
        << " MCTruthStrawHits: " << hits_truth->size()
        << " MCPtr: " << hits_mcptr->size();
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

    size_t nStrawPerEvent = hits->size();
    if (nStrawPerEvent > 0) {
      _hStrawEvt->Fill(nStrawPerEvent);
      _hStrawEvtZ1->Fill(nStrawPerEvent);
      _hStrawEvtZ2->Fill(nStrawPerEvent);
    }
    for (size_t i=0; i<nStrawPerEvent; ++i) {

      // Access data
      StrawHit             const&      hit(hits->at(i));
      StrawHitMCTruth      const&    truth(hits_truth->at(i));
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

      //u coordinate (radial from the center)
      double u = sqrt((xc*xc)+(yc*yc));

      //time of the hit
      double hitTime = hit.time();


      //direction of the straw
      const CLHEP::Hep3Vector stDirection3 = str.getDirection();

      //here the rates
      _hStrawRates[did]->Fill(sid.getStraw(),2*(secid.getSector())+lid.getLayer());
      _hRateUZ->Fill(u,z);
      _hRateU->Fill(u);
      _hRateUT->Fill(u,hitTime);

      // Get MC truth data
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      //Position along the wire using mctruth info
      double vMC = truth.distanceToMid();

      //Position along the wire using dt and propagation velocity (c)
      const double signalVelocity = 299.792458; // mm/ns
      double v = 10e4 * hit.dt()/(2*signalVelocity);

      const CLHEP::Hep3Vector HitPoint = stMidPoint3 + (v/stDirection3.mag())*stDirection3;
      const CLHEP::Hep3Vector MCHitPoint = stMidPoint3 + (vMC/stDirection3.mag())*stDirection3;

      if (fabs(v) > str.getHalfLength()) {
        if (_diagLevel > 0) cout << "Position along the wire bigger than halflength" << endl;
      }

      size_t nHitsPerStraw = mcptr.size();
      _hHitMult->Fill(nHitsPerStraw);

      float tntpArray[56];
      int idx(0);
      tntpArray[idx++] = evt.id().event(); //leaf 1
      tntpArray[idx++] = evt.run(); //leaf 2
      tntpArray[idx++] = hitTime; //leaf 3
      tntpArray[idx++] = hit.dt(); //leaf 4
      tntpArray[idx++] = hitEnergy; //leaf 5
      tntpArray[idx++] = lid.getLayer(); //leaf 6
      tntpArray[idx++] = did; //leaf 7
      tntpArray[idx++] = secid.getSector(); //leaf 8
      tntpArray[idx++] = sid.getStraw(); //leaf 9
      tntpArray[idx++] = MCHitPoint.getX(); //leaf 10
      tntpArray[idx++] = MCHitPoint.getY(); //leaf 11
      tntpArray[idx++] = v; //leaf 12
      tntpArray[idx++] = vMC; //leaf 13
      tntpArray[idx++] = z; //leaf 14

      //Get related G4 hits to identify the track.


      //Map of track id as key, and vector index as value
      map<SimParticleCollection::key_type , size_t > StrawTracksMap;

      //Vectors of pdgId and GenId of the tracks associated to the strawhit
      vector<int>     PdgIdTracks;
      vector<bool>    IsGenerated;
      vector<int>     StepFromEvaVec;
      vector<double>  TracksP;
      vector<int>     TracksStVolumes;

      //List of trackId and energy deposition
      PairList TracksEDep;

      //common index for vectors
      size_t trackIdx(0);

      SimParticleCollection::key_type firstTrackId = SimParticleCollection::key_type(0); 
      bool notFirstTrack = true;
      CLHEP::Hep3Vector const& strDir = str.direction();
      double strRadius = str.getRadius();
      SimParticleCollection::key_type Dau1Idx = SimParticleCollection::key_type(0); 



      for (size_t j = 0; j < mcptr.size(); ++j) {

        StepPointMC const& mchit = *mcptr[j];

        // The simulated particle that made this hit.
        SimParticleCollection::key_type trackId(mchit.trackId());

        //Find in the map if the track is already stored
        map<SimParticleCollection::key_type , size_t >::iterator it;
        it = StrawTracksMap.find(trackId);

        //if the contributing track id does not exist in the map
        //add an element to the map itself, energy to the list and pdgId and genId to the vectors
        if (it==StrawTracksMap.end()) {

	  if (notFirstTrack) {
	    CLHEP::Hep3Vector const& StartPos = simParticles->at(trackId).startPosition();
	    LinePointPCA lppca(stMidPoint3, strDir, StartPos);
	    double insideDistance = lppca.dca();
	    if (insideDistance >= strRadius) {
	      notFirstTrack = true;
	      firstTrackId = trackId;
	    }
	  }
	  
          //insert track id in the trackId vector
          StrawTracksMap.insert(pair<SimParticleCollection::key_type, size_t>(trackId,trackIdx));

          //insert trackId, and energy deposition in the list
          TracksEDep.push_back(pair<SimParticleCollection::key_type, double>(trackId,mchit.eDep()));

          if ( haveSimPart ){
            SimParticle const& sim = simParticles->at(trackId);

            // PDG Particle Id of the sim particle that made this hit.
            PdgIdTracks.push_back(sim.pdgId());
            IsGenerated.push_back(sim.fromGenerator());
	    TracksP.push_back(sim.startMomentum().vect().mag());
	    TracksStVolumes.push_back(sim.startVolumeIndex());

	    int nsfe = 0;

	    if (sim.fromGenerator()) {
	      StepFromEvaVec.push_back(nsfe);
	    } else {
	      bool findEva = false;
	      SimParticle& baby = const_cast<SimParticle&>(sim);
	      while (!findEva) {
		SimParticle & mommy = const_cast<SimParticle&>(*baby.parent());
		nsfe++;
		if (mommy.fromGenerator() || !baby.hasParent()) {
		  if (trackId==firstTrackId) {
		    Dau1Idx = baby.id();
		  }
		  StepFromEvaVec.push_back(nsfe);
		  findEva = true;
		} else {
		  baby = mommy;
		}
	      }
	    }
          } else if ( !haveSimPart) {
            PdgIdTracks.push_back(0);
            IsGenerated.push_back(false);
	    StepFromEvaVec.push_back(0);
	    TracksP.push_back(0);
	    TracksStVolumes.push_back(0);
          }

          //increment index
          trackIdx++;
        } else if (it != StrawTracksMap.end()) {
          for (PairList::iterator it2 = TracksEDep.begin(); it2 != TracksEDep.end(); ++it2) {
            if (it2->first == trackId) {
              it2->second += mchit.eDep();
            }
          }
        }
      }

      TracksEDep.sort(SortByEnergy);

      int nTrkPerStraw = TracksEDep.size();

      if (nTrkPerStraw > 3) {
        if (_diagLevel > 0) {
          cout << "More than 3 different tracks contribute to the straw:"
               << "\nonly the first three with higher e deposit will be stored" << endl;
        }
      }

      tntpArray[idx++] = nTrkPerStraw; //leaf 15
      int counter = 0;

      for (PairList::reverse_iterator it = TracksEDep.rbegin();
           it != TracksEDep.rend(); ++it) {
        if (counter == 3) break;

        size_t vec_idx = StrawTracksMap[it->first];

        tntpArray[idx++] = it->first.asInt(); //leaf 16 - 23 - 30
        tntpArray[idx++] = PdgIdTracks[vec_idx];//leaf 17 - 24 - 31
        tntpArray[idx++] = it->second;//leaf 18 - 25 - 32
        tntpArray[idx++] = IsGenerated[vec_idx];//leaf 19 - 26 - 33
        tntpArray[idx++] = StepFromEvaVec[vec_idx];//leaf 20 - 27 - 34
	tntpArray[idx++] = TracksP[vec_idx];//leaf 21 - 28 - 35
	tntpArray[idx++] = TracksStVolumes[vec_idx];//leaf 22 - 29 - 36
        counter++;
      }

      //Fill with 0 the rest of the ntupla leaves
      //if there are less than 3 tracks contributing
      //to the straw hit

      for (int add_idx = 0; add_idx < 3 - counter; ++add_idx) {
        tntpArray[idx++] = 0;//leaf 16 - 23 - 30 
        tntpArray[idx++] = 0;//leaf 17 - 24 - 31
        tntpArray[idx++] = 0;//leaf 18 - 25 - 32
        tntpArray[idx++] = 0;//leaf 19 - 26 - 33
        tntpArray[idx++] = 0;//leaf 20 - 27 - 34
        tntpArray[idx++] = 0;//leaf 21 - 28 - 35
        tntpArray[idx++] = 0;//leaf 22 - 29 - 36
      }

      size_t vec_idx = StrawTracksMap[firstTrackId];
      tntpArray[idx++] = firstTrackId.asInt(); //leaf 37
      tntpArray[idx++] = PdgIdTracks[vec_idx];//leaf 38
      tntpArray[idx++] = TracksP[vec_idx];//leaf 39
      tntpArray[idx++] = IsGenerated[vec_idx];//leaf 40
      tntpArray[idx++] = TracksStVolumes[vec_idx];//leaf 41
      tntpArray[idx++] = StepFromEvaVec[vec_idx];//leaf 42

      size_t ngen = genParticles->size();
      if (ngen>1) {
        cout << "The plugin is supposed to analyze single background rates,"
             << "with one generated particle per event"
             << "\nThis event has more than one genparticle. Only the "
             << "first one will be stored" << endl;
      }

      SimParticleCollection::key_type idxInSim = SimParticleCollection::key_type(1);
      SimParticle const& geninSim = simParticles->at(idxInSim);
      if (!geninSim.fromGenerator()) {
	cout << "Watch out. First particle is not from generator. What's happening?" << endl;
      }


      if (ngen > 0) {
        GenParticle const& gen = genParticles->at(0);
        tntpArray[idx++] = gen.generatorId().id();//leaf 43
        tntpArray[idx++] = gen.momentum().vect().mag();//leaf 44
        tntpArray[idx++] = gen.momentum().e();//leaf 45
        tntpArray[idx++] = gen.position().x();//leaf 46
        tntpArray[idx++] = gen.position().y();//leaf 47
        tntpArray[idx++] = gen.position().z();//leaf 48
        tntpArray[idx++] = gen.momentum().cosTheta();//leaf 49
        tntpArray[idx++] = gen.momentum().phi();//leaf 50
        tntpArray[idx++] = gen.time();//leaf 51
      } else if ( ngen == 0 ) {
        tntpArray[idx++] = 0;//leaf 43
        tntpArray[idx++] = 0;//leaf 44
        tntpArray[idx++] = 0;//leaf 45
        tntpArray[idx++] = 0;//leaf 46
        tntpArray[idx++] = 0;//leaf 47
        tntpArray[idx++] = 0;//leaf 48
        tntpArray[idx++] = 0;//leaf 49
        tntpArray[idx++] = 0;//leaf 50
        tntpArray[idx++] = 0;//leaf 51
      }


      if (Dau1Idx != SimParticleCollection::key_type(0)) {
	SimParticle const& Dau1 = simParticles->at(Dau1Idx);
	tntpArray[idx++] = Dau1.pdgId();//leaf 52
	tntpArray[idx++] = Dau1.startMomentum().vect().mag();//leaf 53
	tntpArray[idx++] = Dau1.startVolumeIndex();//leaf 54      
      } else {
	tntpArray[idx++] = 0;//leaf 52
	tntpArray[idx++] = 0;//leaf 53
	tntpArray[idx++] = 0;//leaf 54
      }

      tntpArray[idx++] = driftTime; //leaf 55
      tntpArray[idx++] = driftDistance; //leaf 56

      _tNtup->Fill(tntpArray);

    } //end of Strawhits loop

  } // end of doTracker

  void BkgRates::doITracker(art::Event const& evt, bool skip) {

    if (skip) return;

    const Tracker& tracker = getTrackerOrThrow();
    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in PtrStepPointMCVectorCollection, rather than
    // use producer name directly (_g4ModuleLabel).

    if (!(hits->size() == hits_truth->size() &&
          hits_mcptr->size() == hits->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits: " << hits->size()
        << " MCTruthStrawHits: " << hits_truth->size()
        << " MCPtr: " << hits_mcptr->size();
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

    size_t nStrawPerEvent = hits->size();
    if (nStrawPerEvent > 0) {
      _hStrawEvt->Fill(nStrawPerEvent);
      _hStrawEvtZ1->Fill(nStrawPerEvent);
      _hStrawEvtZ2->Fill(nStrawPerEvent);
    }

    for (size_t i=0; i<nStrawPerEvent; ++i) {

      // Access data
      StrawHit             const&  hit(hits->at(i));
      StrawHitMCTruth      const&  truth(hits_truth->at(i));
      PtrStepPointMCVector const&  mcptr(hits_mcptr->at(i));

      double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      if (hitEnergy < _minimumEnergy) continue;

      //Get hit straw
      StrawIndex si = hit.strawIndex();
      const Straw & str = tracker.getStraw(si);
      const Cell & cell = static_cast<const Cell&>( str );


      // cout << "Getting informations about cells" << endl;

      int sid = cell.Id().getCell();
      int lid = cell.Id().getLayer();
      int did = cell.Id().getLayerId().getSuperLayer();

      const CLHEP::Hep3Vector stMidPoint3 = str.getMidPoint();

      //time of the hit
      double hitTime = hit.time();

      //direction of the straw
      const CLHEP::Hep3Vector stDirection3 = str.getDirection();

      // cout << "Reading MCtruth info" << endl;

      // Get MC truth data
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      //Position along the wire using mctruth info
      double vMC = truth.distanceToMid();

      const CLHEP::Hep3Vector MCHitPoint = stMidPoint3 + (vMC/stDirection3.mag())*stDirection3;

      size_t nHitsPerStraw = mcptr.size();
      _hHitMult->Fill(nHitsPerStraw);

      //  cout << "Filling ntupla" << endl;

      float tntpArray[34];
      int idx(0);
      tntpArray[idx++] = evt.id().event();
      tntpArray[idx++] = evt.run();
      tntpArray[idx++] = hitTime;
      tntpArray[idx++] = hitEnergy;
      tntpArray[idx++] = lid;
      tntpArray[idx++] = did;
      tntpArray[idx++] = sid;
      tntpArray[idx++] = MCHitPoint.getX();
      tntpArray[idx++] = MCHitPoint.getY();
      tntpArray[idx++] = vMC;


      //Get related G4 hits to identify the track.

      // cout << "Tracks info" << endl;


      //Map of track id as key, and vector index as value
      map<SimParticleCollection::key_type , size_t > StrawTracksMap;

      //Vectors of pdgId and GenId of the tracks associated to the strawhit
      vector<int>     PdgIdTracks;
      vector<bool>    IsGenerated;

      //List of trackId and energy deposition
      PairList TracksEDep;

      //common index for vectors
      size_t trackIdx(0);

      for (size_t j = 0; j < mcptr.size(); ++j) {

        StepPointMC const& mchit = *mcptr[j];

        // The simulated particle that made this hit.
        SimParticleCollection::key_type trackId(mchit.trackId());

        //Find in the map if the track is already stored
        map<SimParticleCollection::key_type , size_t >::iterator it;
        it = StrawTracksMap.find(trackId);

        //if the contributing track id does not exist in the map
        //add an element to the map itself, energy to the list and pdgId and genId to the vectors
        if (it==StrawTracksMap.end()) {

          //insert track id in the trackId vector
          StrawTracksMap.insert(pair<SimParticleCollection::key_type, size_t>(trackId,trackIdx));

          //insert trackId, and energy deposition in the list
          TracksEDep.push_back(pair<SimParticleCollection::key_type, double>(trackId,mchit.eDep()));

          if ( haveSimPart ){
            SimParticle const& sim = simParticles->at(trackId);

            // PDG Particle Id of the sim particle that made this hit.
            PdgIdTracks.push_back(sim.pdgId());
            IsGenerated.push_back(sim.fromGenerator());

          } else if ( !haveSimPart) {
            PdgIdTracks.push_back(0);
            IsGenerated.push_back(false);
          }

          //increment index
          trackIdx++;
        } else if (it != StrawTracksMap.end()) {
          for (PairList::iterator it2 = TracksEDep.begin(); it2 != TracksEDep.end(); ++it2) {
            if (it2->first == trackId) {
              it2->second += mchit.eDep();
            }
          }
        }
      }

      TracksEDep.sort(SortByEnergy);

      int nTrkPerStraw = TracksEDep.size();

      if (nTrkPerStraw > 3) {
        if (_diagLevel > 0) {
          cout << "More than 3 different tracks contribute to the straw:"
               << "\nonly the first three with higher e deposit will be stored" << endl;
        }
      }

      tntpArray[idx++] = nTrkPerStraw;
      int counter = 0;

      for (PairList::reverse_iterator it = TracksEDep.rbegin();
           it != TracksEDep.rend(); ++it) {
        if (counter == 3) break;

        size_t vec_idx = StrawTracksMap[it->first];

        tntpArray[idx++] = it->first.asInt();
        tntpArray[idx++] = PdgIdTracks[vec_idx];
        tntpArray[idx++] = it->second;
        tntpArray[idx++] = IsGenerated[vec_idx];
        counter++;
      }


      //Fill with 0 the rest of the ntupla leaves
      //if there are less than 3 tracks contributing
      //to the straw hit

      for (int add_idx = 0; add_idx < 3 - counter; ++add_idx) {
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
      }





      size_t ngen = genParticles->size();
      if (ngen>1) {
        cout << "The plugin is supposed to analyze single background rates,"
             << "with one generated particle per event"
             << "\nThis event has more than one genparticle. Only the "
             << "first one will be stored" << endl;
      }
      if (ngen > 0) {
        GenParticle const& gen = genParticles->at(0);
        tntpArray[idx++] = gen.generatorId().id();
        tntpArray[idx++] = gen.momentum().vect().mag();
        tntpArray[idx++] = gen.momentum().e();
        tntpArray[idx++] = gen.position().x();
        tntpArray[idx++] = gen.position().y();
        tntpArray[idx++] = gen.position().z();
        tntpArray[idx++] = gen.momentum().cosTheta();
        tntpArray[idx++] = gen.momentum().phi();
        tntpArray[idx++] = gen.time();
      } else if ( ngen == 0 ) {
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
        tntpArray[idx++] = 0;
      }

      tntpArray[idx++] = driftTime;
      tntpArray[idx++] = driftDistance;

      _tNtup->Fill(tntpArray);

    } //end of Cellhits loop

  } // end of doITracker

  void BkgRates::doCalorimeter(art::Event const& evt, bool skip) {

    if (skip) return;

    const double CrDensity = 8.28; //from G4. It is in g/cm^3

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    double CrSize = cg->crystalHalfSize();
    double CrLength = cg->crystalHalfLength();
    double CrVolumeCm3 = (CrSize*CrSize*CrLength)/1000; //factor 1000 because they are in mm

    double CrMassKg = CrDensity*CrVolumeCm3/1000;

    //    cout << CrSize << '\t' << CrLength << '\t' << CrVolumeCm3
    //     << '\t' << CrDensity << '\t' << CrMassKg << endl;

    // Get handles to calorimeter collections
    art::Handle<CaloHitCollection> caloHits;
    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    //    art::Handle<StepPointMCCollection> steps;

    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);
    //    evt.getByLabel(_g4ModuleLabel,"calorimeter",steps);
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
      _hCryEvt->Fill(caloCrystalHits->size());
      _hCryEvtZ1->Fill(caloCrystalHits->size());
      _hCryEvtZ2->Fill(caloCrystalHits->size());

      for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {

        CaloCrystalHit const & hit = (*caloCrystalHits).at(i);

        std::vector<art::Ptr<CaloHit> > const & ROIds  = hit.readouts();

	//	cout << "Event " << evt.id().event() << ". In the caloCrystalHits there are " << ROIds.size() << " RO associated" << endl;

        if (ROIds.size() < 1) {
          //          cout << " Event n. " << evt.id().event()
          //   << " \t got crystal hits but not ROhits"
          //     << '\t' <<  caloCrystalHits->size()
          //   << '\t' << ROIds.size() << endl;
          continue;
        }

        double EfromShower = 0.;
        double EfromOutside1 = 0.;
        double EfromOutside2 = 0.;
        double EfromOutside3 = 0.;
        double EfromOutsideAll = 0.;
        double EfromOtherVane = 0.;
        int OutsideTrkPdgId1 = 0;
        int OutsideTrkPdgId2 = 0;
        int OutsideTrkPdgId3 = 0;
        int nOutsideTrk = 0;
        double GeneratedEDep = 0.;

        if (hit.energyDep() < _minimumEnergy) continue;

        bool readCryOnce(false);
        float cntpArray[36];
        int idx(0);

        //List of trackId and energy deposition
        PairList OutsideEDep;

	double firstHitTime = 100000;
	CLHEP::Hep3Vector firstHitPos(0,0,0);

	CLHEP::Hep3Vector cryFrame(0,0,0);


        for (size_t it = 0;
             it < ROIds.size() ; ++it ) {

          size_t collectionPosition = ROIds.at(it).key();
          CaloHit const & thehit = *ROIds.at(it);

	  //	  cout << "ROID n. " << it << ": informations " << (readCryOnce?"not":"") << " stored.\n"
	  //       << "Crystal : " << cg->getCrystalByRO(thehit.id())   
          //     << "\nEnergy: " << hit.energyDep() << endl;
 

          if (!readCryOnce) {
            CLHEP::Hep3Vector cryCenter =  cg->getCrystalOriginByRO(thehit.id());
            int vane = cg->getVaneByRO(thehit.id());
            int Zcry = cg->getCrystalZByRO(thehit.id());
            int Rcry = cg->getCrystalRByRO(thehit.id());
            _hCrystalRates[vane]->Fill(Zcry,Rcry);
            cntpArray[idx++] = evt.id().event();
            cntpArray[idx++] = evt.run();
            cntpArray[idx++] = hit.time();
            cntpArray[idx++] = hit.energyDep();
            double CrEjoule = hit.energyDep()*CLHEP::joule/CLHEP::megaelectronvolt;
            double dose = CrEjoule/CrMassKg;
            cntpArray[idx++] = dose;
            cntpArray[idx++] = cg->getCrystalByRO(thehit.id());
            cntpArray[idx++] = vane;
            cntpArray[idx++] = cryCenter.getX() + 3904.;  //value used to shift in tracker coordinate system
            cntpArray[idx++] = cryCenter.getY();
            cntpArray[idx++] = cryCenter.getZ() - 10200;  //value used to shift in tracker coordinate system
	  	  
	    
	    PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
	    size_t nHitsPerCrystal = mcptr.size();
	    _hCaloHitMult->Fill(nHitsPerCrystal);
	    
	    //	    cout << "In the RO there are " << nHitsPerCrystal << " hits. List index is " << collectionPosition << endl;
	    
	    for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {
	      
	      StepPointMC const& mchit = *mcptr[j2];
	      // The simulated particle that made this hit.
	      SimParticleCollection::key_type trackId(mchit.trackId());
	      {
		SimParticle const& sim = simParticles->at(trackId);
		if (sim.stoppingCode() == ProcessCode::mu2eMaxSteps) {
		  cout << "Track " << sim.pdgId() << "  dies at "
		       << sim.endGlobalTime() << endl;
		}
	      }
	      
	      CaloManager->setTrackAndRO(evt, _g4ModuleLabel, trackId, thehit.id() );
	      
	      //	      cout << "Original Vane: " << CaloManager->localVane()
	      //		   << "\nStarting Vane: " << CaloManager->startingVane() << endl;
	      
	      if (CaloManager->localVane() == CaloManager->startingVane()) {
		EfromShower += mchit.eDep();
		//		cout << "From shower we have " << mchit.eDep() << " and globally " << EfromShower << endl;
	      }
	      
	      if (CaloManager->localVane() != CaloManager->startingVane() &&
		  CaloManager->startingVane() != -1) {
		EfromOtherVane += mchit.eDep();
		//		cout << "From other vane we have " << mchit.eDep() << " and globally " << EfromOtherVane << endl;
	      }
	      
	      if (CaloManager->fromOutside()) {
		if (!(CaloManager->generated())) {
		  PairListAdd(OutsideEDep, trackId, mchit.eDep());
		  EfromOutsideAll += mchit.eDep();
		  //		  cout << "From outside we have " << mchit.eDep() << " and globally " << EfromOutsideAll << endl;
		}
		
		if (CaloManager->generated()) {
		  GeneratedEDep +=  mchit.eDep();
		  if (mchit.time()<firstHitTime) {
		    firstHitTime = mchit.time();
		    firstHitPos = mchit.position();
		    //		    cout << "before " << firstHitPos << endl;
		    cryFrame = cg->toCrystalFrame(thehit.id(), firstHitPos);
		    //		    cout << "after " << firstHitPos << endl;
		  }
		  //		  cout << "From generated we have " << mchit.eDep() << " and globally " << GeneratedEDep << endl;
		  //		  cout << "Time of this hit is " << mchit.time() << " and position is " << mchit.position() << endl;
		}
	      }
	    }
	  }
	  
	  readCryOnce = true;
	  
	}
	
	OutsideEDep.sort(SortByEnergy);

        nOutsideTrk = OutsideEDep.size();
	
        if (nOutsideTrk > 3) {
          if (_diagLevel > 0) {
            cout << "More than 3 different tracks from outside the calorimeter contribute to the crystal:"
                 << "\nonly the first three with higher e deposit will be stored" << endl;
          }
        }

        int counter = 0;

        PairList::reverse_iterator it = OutsideEDep.rbegin();

        if (nOutsideTrk>0) {

          for (int itidx=0; itidx<nOutsideTrk; ++itidx) {
            if (counter == 3) break;
            if (itidx == 0) {
              EfromOutside1 = it->second;
              SimParticle const& sim = simParticles->at(it->first);
              OutsideTrkPdgId1 = sim.pdgId();
              it++;
            }
            if (itidx == 1) {
              EfromOutside2 = it->second;
              SimParticle const& sim = simParticles->at(it->first);
              OutsideTrkPdgId2 = sim.pdgId();
              it++;
            }
            if (itidx == 2) {
              EfromOutside3 = it->second;
              SimParticle const& sim = simParticles->at(it->first);
              OutsideTrkPdgId3 = sim.pdgId();
              it++;
            }
          }

        }

        cntpArray[idx++] = EfromShower;
        cntpArray[idx++] = EfromOtherVane;
        cntpArray[idx++] = nOutsideTrk;
        cntpArray[idx++] = EfromOutside1;
        cntpArray[idx++] = OutsideTrkPdgId1;
        cntpArray[idx++] = EfromOutside2;
        cntpArray[idx++] = OutsideTrkPdgId2;
        cntpArray[idx++] = EfromOutside3;
        cntpArray[idx++] = OutsideTrkPdgId3;
	cntpArray[idx++] = EfromOutsideAll;
        cntpArray[idx++] = GeneratedEDep;
	//	cout << "The one I store has the following position " << firstHitPos << endl;
	cntpArray[idx++] = firstHitPos.x();
	cntpArray[idx++] = firstHitPos.y();
	cntpArray[idx++] = firstHitPos.z();
	cntpArray[idx++] = cryFrame.x();
	cntpArray[idx++] = cryFrame.y();
	cntpArray[idx++] = cryFrame.z();

        size_t ngen = genParticles->size();
        if (ngen>1) {
          cout << "The plugin is supposed to analyze single background rates,"
               << "with one generated particle per event"
               << "\nThis event has more than one genparticle. Only the "
               << "first one will be stored" << endl;
        }
        if (ngen > 0) {
          GenParticle const& gen = genParticles->at(0);
          cntpArray[idx++] = gen.generatorId().id();
          cntpArray[idx++] = gen.momentum().vect().mag();
          cntpArray[idx++] = gen.momentum().e();
          cntpArray[idx++] = gen.position().x();
          cntpArray[idx++] = gen.position().y();
          cntpArray[idx++] = gen.position().z();
          cntpArray[idx++] = gen.momentum().cosTheta();
          cntpArray[idx++] = gen.momentum().phi();
          cntpArray[idx++] = gen.time();
        } else if ( ngen == 0 ) {
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
          cntpArray[idx++] = 0.;
        }
        if (hit.energyDep()>0) {
          _cNtup->Fill(cntpArray);
        }
      }
    }
  } // end of doCalorimeter

  void BkgRates::doStoppingTarget(const art::Event& event) {

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

      if ( volInfo.name() == "TargetFoil_" ) {

        if( stoppedtracks.insert(trackId).second ) {
          float tgtntpArray[10];
          int idx(0);
          tgtntpArray[idx++] = event.id().event();
          tgtntpArray[idx++] = event.run();
          tgtntpArray[idx++] = sim->endGlobalTime();
          tgtntpArray[idx++] = sim->endPosition().x();
          tgtntpArray[idx++] = sim->endPosition().y();
          tgtntpArray[idx++] = sim->endPosition().z();
          tgtntpArray[idx++] = sim->fromGenerator();
          if (sim->fromGenerator()) generatedStopped = true;
          tgtntpArray[idx++] = sim->pdgId();
          tgtntpArray[idx++] = trackId.asInt();
          tgtntpArray[idx++] = volInfo.copyNo();

          _tgtNtup->Fill(tgtntpArray);

        }
      }
    }
    _skipEvent = generatedStopped;
  }  // end doStoppingTarget
}

using mu2e::BkgRates;
DEFINE_ART_MODULE(BkgRates);

