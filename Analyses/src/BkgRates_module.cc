// C++ includes
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "Analyses/inc/MCCaloUtilities.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloCrystalOnlyHitCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/GenParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

// Other external includes
#include "CLHEP/Units/PhysicalConstants.h"

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
      _nZCryPerVane(50),
      _nRCryPerVane(10)
    {
    }
    virtual ~BkgRates() {
    }
    virtual void beginJob();

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

    art::ServiceHandle<GeometryService> geom;

    if (ncalls == 1) {

      // cout << "This should be done only in the first event" << endl;


      art::ServiceHandle<art::TFileService> tfs;

      if (geom->hasElement<TTracker>()) {
        _tNtup        = tfs->make<TNtuple>( "StrawHits", "Straw Ntuple",
                                            "evt:run:time:dt:eDep:lay:dev:sec:strawId:MChitX:MChitY:v:vMC:z:nTrk:t1trkId:t1pdgId:t1en:t1isGen:t2trkId:t2pdgId:t2en:t2isGen:t3trkId:t3pdgId:t3en:t3isGen:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime:driftTime:driftDist" );
      } else if(geom->hasElement<ITracker>()) {
        _tNtup        = tfs->make<TNtuple>( "CellHits", "Cell Ntuple",
                                            "evt:run:time:eDep:lay:superlay:cellId:MChitX:MChitY:wireZMC:nTrk:t1trkId:t1pdgId:t1en:t1isGen:t2trkId:t2pdgId:t2en:t2isGen:t3trkId:t3pdgId:t3en:t3isGen:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime:driftTime:driftDist" );
      }
      _cNtup        = tfs->make<TNtuple>( "CaloHits", "Calo Ntuple",
                                          "evt:run:crTime:crE:crRad:crId:crVane:crX:crY:crZ:ESwr:EOutVane:NtrkOutside:OutsideE1:OutsidePdg1:OutsideE2:OutsidePdg2:OutsideE3:OutsidePdg3:EGen:genId:genP:genE:genX:genY:genZ:genCosTh:genPhi:genTime" );
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
    art::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than
    // use producer name directly (_g4ModuleLabel).

    art::Handle<StepPointMCCollection> mchitsHandle;
    evt.getByLabel(_g4ModuleLabel,_trackerStepPoints,mchitsHandle);
    StepPointMCCollection const* mchits = mchitsHandle.product();

    /*
    //Get and printout all straw (just to check)
    const std::deque<Straw>& allstraws = tracker.getAllStraws();
    for (size_t i = 0; i< allstraws.size(); ++i) {
    Straw str = allstraws[i];
    cout << "Index " << str.Index()
    << "\t" << str.Id() << endl;
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
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));

      double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      if (hitEnergy < _minimumEnergy) continue;

      //Get hit straw
      StrawIndex si = hit.strawIndex();
      Straw str = tracker.getStraw(si);
      StrawId sid = str.Id();
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

      float tntpArray[38];
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
      tntpArray[idx++] = MCHitPoint.getX();
      tntpArray[idx++] = MCHitPoint.getY();
      tntpArray[idx++] = v;
      tntpArray[idx++] = vMC;
      tntpArray[idx++] = z;

      //Get related G4 hits to identify the track.


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

        StepPointMC const& mchit = (*mchits)[mcptr[j].index];

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
    art::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than
    // use producer name directly (_g4ModuleLabel).

    art::Handle<StepPointMCCollection> mchitsHandle;
    evt.getByLabel(_g4ModuleLabel,_trackerStepPoints,mchitsHandle);
    StepPointMCCollection const* mchits = mchitsHandle.product();

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
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));

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

        StepPointMC const& mchit = (*mchits)[mcptr[j].index];

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
    art::Handle<CaloHitMCTruthCollection> caloMC; //unused
    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<DPIndexVectorCollection> mcptrHandle;
    art::Handle<StepPointMCCollection> steps;

    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);
    evt.getByLabel(_g4ModuleLabel,"calorimeter",steps);
    evt.getByLabel(_caloReadoutModuleLabel, caloHits);
    evt.getByLabel(_caloReadoutModuleLabel, caloMC);
    evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();
    StepPointMCCollection const* mchits = steps.product();
    if (!( caloHits.isValid() && caloMC.isValid())) {
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

    map<unsigned, size_t> CHMap;

    for (size_t j=0; j<caloHits->size(); ++j) {

      CaloHit const & ahit = (*caloHits).at(j);
      CHMap.insert(pair<unsigned, size_t>(ahit.id(), j));

    }

    if (caloCrystalHits->size()>0) {
      _hCryEvt->Fill(caloCrystalHits->size());
      _hCryEvtZ1->Fill(caloCrystalHits->size());
      _hCryEvtZ2->Fill(caloCrystalHits->size());

      for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {

        CaloCrystalHit const & hit = (*caloCrystalHits).at(i);

        DPIndexVector const & ROIds  = hit.roIds();

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
        double EfromOtherVane = 0.;
        int OutsideTrkPdgId1 = 0;
        int OutsideTrkPdgId2 = 0;
        int OutsideTrkPdgId3 = 0;
        int nOutsideTrk = 0;
        double GeneratedEDep = 0.;

        if (hit.energyDep() < _minimumEnergy) continue;

        bool readCryOnce(false);
        float cntpArray[29];
        int idx(0);

        //List of trackId and energy deposition
        PairList OutsideEDep;

        for (size_t it = 0;
             it < ROIds.size() ; ++it ) {

          size_t CollectionPosition = CHMap[ROIds.at(it).index];

          CaloHit const & thehit = (*caloHits).at(CollectionPosition);

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


            DPIndexVector const & mcptr(hits_mcptr->at(CollectionPosition));
            size_t nHitsPerCrystal = mcptr.size();
            _hCaloHitMult->Fill(nHitsPerCrystal);

            for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

              StepPointMC const& mchit = (*mchits)[mcptr[j2].index];
              // The simulated particle that made this hit.
              SimParticleCollection::key_type trackId(mchit.trackId());

              CaloManager->setTrackAndRO(evt, _g4ModuleLabel, trackId, ROIds.at(it).index );

              //cout << "Original Vane: " << CaloManager->localVane()
              //     << "\nStarting Vane: " << CaloManager->startingVane() << endl;
              if (CaloManager->localVane() == CaloManager->startingVane()) {
                EfromShower += mchit.eDep();
              }
              if (CaloManager->localVane() != CaloManager->startingVane() &&
                  CaloManager->startingVane() != -1) {
                EfromOtherVane += mchit.eDep();
              }

              if (CaloManager->fromOutside()) {
                if (!CaloManager->generated()) {
                  PairListAdd(OutsideEDep, trackId, mchit.eDep());
                }

                if (CaloManager->generated()) {
                  GeneratedEDep +=  mchit.eDep();
                }
              }
            }

            readCryOnce = true;

          }
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
        cntpArray[idx++] = GeneratedEDep;

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

      SimParticle const* sim = simParticles->findOrNull(trackId);
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

