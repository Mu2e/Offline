// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <deque>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Provenance/interface/Provenance.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloCrystalHitMCTruthCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

using namespace std;

namespace mu2e {


  class BkgRates : public edm::EDAnalyzer {
  public:
    explicit BkgRates(edm::ParameterSet const& pset):
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel")),
      _nAnalyzed(0),
      _hHitMult(0),
      _hStrawEvt(0),
      _hRateUZ(0),
      _hRateU(0),
      _hRateUT(0),
      _hRateMaxU(0),
      _tNtup(0),
      _cNtup(0),
      _nGenParticles(0),
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
      //cout << "NgenParticles " << _nGenParticles << endl;

      for (int i=0; i<_nDevices; ++i) {
        _hStrawRates[i]->Scale(1/(double)_nGenParticles);
      }

      for (int i=0; i<_nVanes; ++i) {
        _hCrystalRates[i]->Scale(1/(double)_nGenParticles);
      }

      _hRateUZ->Scale(1/(double)_nGenParticles);
      _hRateU->Scale(1/(double)_nGenParticles);
      _hRateUT->Scale(1/(double)_nGenParticles);

    }
    
    virtual void beginJob(edm::EventSetup const&);

    void analyze(edm::Event const& e, edm::EventSetup const&);

  private:

    void doTracker(edm::Event const& evt);

    void doCalorimeter(edm::Event const& evt);

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    //number of analyzed events
    int _nAnalyzed;


    TH1F* _hHitMult;
    TH1F* _hStrawEvt;
    TH2F* _hRateUZ;
    TH1F* _hRateU;
    TH2F* _hRateUT;
    TH1F* _hRateMaxU;
    
    TNtuple* _tNtup;
    TNtuple* _cNtup;

    int _nGenParticles;

    const int _nDevices, _nSectors, _nLayers, _nStrawsPerLay;
    const int _nVanes, _nZCryPerVane, _nRCryPerVane;

    vector<TH2F*> _hStrawRates; 
    vector<TH2F*> _hCrystalRates;

  };

  void BkgRates::beginJob(edm::EventSetup const& ) {

    edm::Service<edm::TFileService> tfs;

    _hHitMult     = tfs->make<TH1F>( "hHitMult",    "Multiplicity of g4 hit per Straw ", 100,  0.,  100. );
    _hStrawEvt    = tfs->make<TH1F>( "hStrawEvt",   "Multiplicity of straws per event ", 200,  0.,  200. );
    _hRateUZ      = tfs->make<TH2F>( "hRateUZ",     "Straw hit in u and z coordinates",  80,  380.,  690., 40, -1550., 1550.);
    _hRateU       = tfs->make<TH1F>( "hRateU",      "Straw hit in u coordinate",         100,  380.,  690. ); 
    _hRateUT      = tfs->make<TH2F>( "hRateUT",     "Straw hit in u coordinate and time",80,  380.,  690., 40, 700., 1900.);
    //    _hRateMaxU    = tfs->make<TH1F>( "hRateMaxU",      ask Asset

    _tNtup        = tfs->make<TNtuple>( "StrawHits", "Straw Ntuple",
                                        "evt:time:dt:eDep:lay:dev:sec:strawId:hLeng:dirX:dirY:dirZ:x:y:u:v:z:wireX:wireY:pdgId:genId" );
    _cNtup        = tfs->make<TNtuple>( "CaloHits", "Calo Ntuple",
                                        "evt:crTime:crE:crId:hTime:hE:hId:crX:crY:crZ" );

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
  
  void BkgRates::analyze(edm::Event const& evt, edm::EventSetup const&) {

    ++_nAnalyzed;

    //*****test code******
        static int ncalls(0);
    ++ncalls;
    if (ncalls == 1) {
      // cout << "This should be done only in the first event" << endl;
    }

    doTracker(evt);

    doCalorimeter(evt);

  } // end of analyze

  void BkgRates::doTracker(edm::Event const& evt) {

    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    edm::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.

    edm::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs

    edm::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than 
    // use producer name directly ("g4run"). 

    edm::Handle<StepPointMCCollection> mchitsHandle;
    evt.getByLabel("g4run",_trackerStepPoints,mchitsHandle);
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
      throw cms::Exception("RANGE")
        << "Strawhits: " << hits->size() 
        << " MCTruthStrawHits: " << hits_truth->size() 
        << " MCPtr: " << hits_mcptr->size(); 
    }

    // Get handles to the generated and simulated particles.
    edm::Handle<ToyGenParticleCollection> genParticles;
    evt.getByType(genParticles);

    _nGenParticles += genParticles->size();

    edm::Handle<SimParticleCollection> simParticles;
    evt.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    edm::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    size_t nStrawPerEvent = hits->size();
    _hStrawEvt->Fill(nStrawPerEvent);

    for (size_t i=0; i<nStrawPerEvent; ++i) {
      
      // Access data
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));

      //Get hit straw
      StrawIndex si = hit.strawIndex();
      Straw str = tracker.getStraw(si);
      StrawId sid = str.Id();
      LayerId lid = sid.getLayerId();
      DeviceId did = sid.getDeviceId();
      SectorId secid = sid.getSectorId();
      
      const CLHEP::Hep3Vector stMidPoint3 = str.getMidPoint();
      const CLHEP::Hep3Vector stDirection3 = str.getDirection();
      double xc = stMidPoint3.getX();
      double yc = stMidPoint3.getY();
      double u = sqrt((xc*xc)+(yc*yc));
      double z = stMidPoint3.getZ();

      //here the rates
      _hStrawRates[did]->Fill(sid.getStraw(),2*(secid.getSector())+lid.getLayer());

      // Use MC truth data
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();
      double distanceToMid = truth.distanceToMid();

      size_t nHitsPerStraw = mcptr.size();
      _hHitMult->Fill(nHitsPerStraw);

      //Get related G4 hits
      for (size_t j = 0; j < mcptr.size(); ++j) {

        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
        
        //Get v coordinate (radial to the center of the panel)
        //Check if x,y origin are in the center of the panel
        
        const CLHEP::Hep3Vector& hitPos3 = mchit.position();

        double v = stDirection3.dot(hitPos3-stMidPoint3);

        const CLHEP::Hep3Vector& point  = stMidPoint3 + v*stDirection3.unit();

        // The simulated particle that made this hit.
        SimParticleCollection::key_type trackId(mchit.trackId());
        
        // Default values for these, in case information is not available.
        int pdgId(0);
        GenId genId;

        if ( haveSimPart ){
          SimParticle const& sim = simParticles->at(trackId);

          // PDG Particle Id of the sim particle that made this hit.
          pdgId = sim.pdgId();
      
          // If this is a generated particle, which generator did it come from?
          // This default constructs to "unknown".
          if ( sim.fromGenerator() ){
            ToyGenParticle const& gen = genParticles->at(sim.generatorIndex());
            genId = gen._generatorId;
          } 
        } 

        float hitTime = hit.time();

        _hRateUZ->Fill(u,z);
        _hRateU->Fill(u);
        _hRateUT->Fill(u,hitTime);

        float tntpArray[23];
        int idx(0);
        tntpArray[idx++] = evt.id().event();
        tntpArray[idx++] = hitTime;
        tntpArray[idx++] = hit.dt();
        tntpArray[idx++] = hit.energyDep();
        tntpArray[idx++] = lid.getLayer();
        tntpArray[idx++] = did;
        tntpArray[idx++] = secid.getSector();
        tntpArray[idx++] = sid.getStraw();
        tntpArray[idx++] = str.getHalfLength();
        tntpArray[idx++] = stDirection3.getX();
        tntpArray[idx++] = stDirection3.getY();
        tntpArray[idx++] = stDirection3.getZ();
        tntpArray[idx++] = xc;
        tntpArray[idx++] = yc;
        tntpArray[idx++] = u;
        tntpArray[idx++] = v;
        tntpArray[idx++] = z;
        tntpArray[idx++] = point.getX();
        tntpArray[idx++] = point.getY();
        tntpArray[idx++] = pdgId;
        tntpArray[idx++] = genId.Id();
        tntpArray[idx++] = driftTime;
        tntpArray[idx++] = driftDistance;
        tntpArray[idx++] = distanceToMid;
      


        _tNtup->Fill(tntpArray);

      } // end of steppointhits loop

    } //end of Strawhits loop

  } // end of doTracker

  void BkgRates::doCalorimeter(edm::Event const& evt) {

    //Get handle to the calorimeter
    edm::Service<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;


    //Test printout of Calorimeter geometry and classification
    /*   
         cout << "Calorimeter geometry: " << endl;
         //    for (int nv=0; nv<cg->nVane(); ++nv) {
         //  for (int nc=0; nc<cg->nCrystalPerVane(); ++nc) {
         int iCry = -1;
         for (size_t nRO=0; nRO<cg->nRO(); ++nRO) {
         int tiCry = cg->getCrystalByRO(nRO);
         if (tiCry != iCry) {
         cout << "\nVane n. " << cg->getVaneByRO(nRO)
         << "\tCrystal (within a vane) n. " << cg->getCrystalVaneByRO(nRO)
         << "\tCrystal n. " << cg->getCrystalByRO(nRO) 
         << "\nPosition " << cg->getCrystalOriginByRO(nRO)
         << "\t RO n. " << flush;
         iCry = tiCry;
         }
         cout << '\t' << nRO << flush;  
         }
    */

    // Get handles to calorimeter collections
    edm::Handle<CaloHitCollection> caloHits;
    edm::Handle<CaloHitMCTruthCollection> caloMC; //unused
    edm::Handle<CaloCrystalHitCollection>  caloCrystalHits;

    evt.getByType(caloHits);
    evt.getByType(caloMC);
    evt.getByType(caloCrystalHits);

    if (!( caloHits.isValid() && caloMC.isValid())) {
      return;
    }

    if (!caloCrystalHits.isValid()) {
      cout << "NO CaloCrystalHits" << endl;
      return;
    }

    typedef multimap<int,size_t> crystalsHitsMultiMap;
    crystalsHitsMultiMap crystalsHits;

    for ( size_t i=0; i<caloHits->size(); ++i ) {

      int cid = cg->getCrystalByRO(caloHits->at(i).id());
      crystalsHits.insert(pair<int,size_t>(cid,i));
      
      /*      cout << "inserting info on map:\n" 
              << "CryId " << cid << '\t' << "hit index " << caloHits->at(i).id() << '\n'
              << "energy " << hit.energyDep() << '\n'
              << "time " << hit.time() << endl;
      */
      
    }

    for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {

      CaloCrystalHit const & hit = (*caloCrystalHits).at(i);
      int cid = hit.id();

      CLHEP::Hep3Vector cryCenter(0,0,0);

      pair<multimap<int,size_t>::iterator,multimap<int,size_t>::iterator> itRange; 
      itRange = crystalsHits.equal_range(cid);

      multimap<int,size_t>::iterator it;

      bool readCryOnce(false);




      for (it = itRange.first; it!= itRange.second; ++it) {

        CaloHit const & clHit = caloHits->at(it->second);  
        
        if (!readCryOnce) {
          cryCenter =  cg->getCrystalOriginByRO(clHit.id());
          int vane = cg->getVaneByRO(clHit.id());
          int Zcry = cg->getCrystalZByRO(clHit.id());
          int Rcry = cg->getCrystalRByRO(clHit.id());
          _hCrystalRates[vane]->Fill(Zcry,Rcry);
          readCryOnce = true;
        }
        
        float cntpArray[10];
        int idx(0);

        cntpArray[idx++] = evt.id().event();
        cntpArray[idx++] = hit.time();
        cntpArray[idx++] = hit.energyDep();
        cntpArray[idx++] = cid;
        cntpArray[idx++] = clHit.time();
        cntpArray[idx++] = clHit.energyDep();
        cntpArray[idx++] = clHit.id();
        cntpArray[idx++] = cryCenter.getX() + 3904.;  //value used to shift in tracker coordinate system
        cntpArray[idx++] = cryCenter.getY();
        cntpArray[idx++] = cryCenter.getZ() - 10200;  //value used to shift in tracker coordinate system

        /*        cout << "form ntupla :\n"
                  << "cry id " << cid << '\n' 
                  << "cry e " << hit.energyDep() << '\n'
                  << "cry time " << hit.time() << '\n'
                  << "caloHit id " << clHit.id() << '\n'
                  << "caloHit e " << clHit.energyDep() << '\n'
                  << "caloHit time " << clHit.time() << endl;
        */

        _cNtup->Fill(cntpArray);        
        
      } // end of loop on hits coming from the same crystal
      
    } // end of loop on crystal hits
    
    
    // Unused caloCrystalHitMCTruths
    edm::Handle<CaloCrystalHitMCTruthCollection>  caloCrystalHitMCTruths;
    evt.getByType(caloCrystalHitMCTruths);
    if (!caloCrystalHitMCTruths.isValid()) {
      cout << ": NO CaloCrystalHitMCTruths" << endl;
      return;
    }
    
    
  } // end of doCalorimeter

}

using mu2e::BkgRates;
DEFINE_FWK_MODULE(BkgRates);
 
