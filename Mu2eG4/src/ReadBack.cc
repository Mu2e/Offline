//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: ReadBack.cc,v 1.44 2011/05/18 21:14:30 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 21:14:30 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes.
#include "Mu2eG4/src/ReadBack.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloCrystalOnlyHitCollection.hh"
#include "ToyDP/inc/DPIndexVector.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StatusG4.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/unknownPDGIdName.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarIndex.hh"

#include "G4Helper/inc/G4Helper.hh"


// Root includes.
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TGraph.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  ReadBack::ReadBack(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _caloCrystalHitsMaker(pset.get<string>("caloCrystalHitsMaker","CaloCrystalHitsMaker")),
    _targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget")),
    _crvStepPoints(pset.get<string>("CRVStepPoints","CRV")),
    _minimumEnergy(pset.get<double>("minimumEnergy")),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _xyHitsMax(pset.get<int>("xyHitsMax",10000)),

    // Histograms
    _nAnalyzed(0),
    _hRadius(0),
    _hTime(0),
    _hMultiplicity(0),
    _hDriftDist(0),
    _hDriftDistW(0),
    _hxHit(0),
    _hyHit(0),
    _hzHit(0),
    _hHitNeighbours(0),
    _hCheckPointRadius(0),
    _hCheckPointRadiusW(0),
    _hCheckPointWireZ(0),
    _hMomentumG4(0),
    _hStepLength(0),
    _hEdep(0),
    _hEdepMC(0),
    _hNcrystal(0),
    _hNcrstep(0),
    _hNrostep(0),
    _hEdepROMC(0),
    _hRCEdep(0),
    _hRCTime(0),
    _hRCNCrystals(0),
    _hRCEdepMC(0),
    _hRCTimeMC(0),
    _hRCNCrystalsMC(0),
    _hTargetEdep(0),
    _hTargetPathLength(0),
    _hTargetNfoils(0),
    _hTargetNfoils2D(0),
    _ntup(0),
    _xyHits(0),
    _xyHitCount(0),
    // CRV
    _hCRVMultiplicity(0),
    _ntupCRV(0),
    //
    // Remaining member data
    _nBadG4Status(0){
  }

  void ReadBack::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius",       "Radius of Hits;(mm)",     100,  0., 1000. );
    _hEnergyDep    = tfs->make<TH1F>( "hEnergyDep",    "Energy Deposited;(keV)",  100,  0.,   10. );
    _hTime         = tfs->make<TH1F>( "hTime",         "Pulse Height;(ns)",       100,  0., 2000. );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event",          100,  0.,  300. );
    _hDriftDist    = tfs->make<TH1F>( "hDriftDist", "Crude Drift Distance;(mm)",  100,  0.,   3.  );
    _hDriftDistW   = tfs->make<TH1F>( "hDriftDistW", "Normalized Crude Drift Distance", 150,  0., 1.5 );

    _hxHit         = tfs->make<TH1F>( "hxHit",  "X of Hit;(mm)",                  100, -1000., 1000. );
    _hyHit         = tfs->make<TH1F>( "hyHit",  "Y of Hit;(mm)",                  100, -1000., 1000. );
    _hzHit         = tfs->make<TH1F>( "hzHit",  "Z of Hit;(mm)",                  100, -1400., 1400. );

    _hHitNeighbours    = tfs->make<TH1F>( "hHitNeighbours",  "Number of hit neighbours",
                                          10, 0., 10. );

    _hCheckPointRadius  = tfs->make<TH1F>( "hCheckPointRadius", "Radius of Reference point; (mm)",
                                           100, 2.25, 2.75 );
    _hCheckPointRadiusW = tfs->make<TH1F>( "hCheckPointRadiusW","Normalized Radius of Reference point",
                                           200, 0., 2. );
    _hCheckPointWireZ   = tfs->make<TH1F>( "hCheckPointWireZ",  "Normalized Wire Z of Reference point",
                                           300, -1.5, 1.5 );

    _hMomentumG4 = tfs->make<TH1F>( "hMomentumG4",  "Mommenta of particles created inside G4; (MeV)",
                                    100, 0., 100. );
    _hStepLength = tfs->make<TH1F>( "hStepLength",  "G4 Step Length in Sensitive Detector; (mm)",
                                    100, 0., 10. );

    _hEdep     = tfs->make<TH1F>( "hEdep",     "Total energy deposition in calorimeter", 2400, 0., 2400. );
    _hEdepMC   = tfs->make<TH1F>( "hEdepMC",   "True energy deposition in calorimeter",  240,  0., 240. );
    _hNcrystal = tfs->make<TH1F>( "hNcrystal", "Total energy deposition in calorimeter",   50, 0.,   50. );
    _hNcrstep  = tfs->make<TH1F>( "hNcrstep", "Number of G4 steps in crystal, per APD", 100, 0., 500. );
    _hNrostep  = tfs->make<TH1F>( "hNrostep", "Number of G4 steps in APD, per APD", 10, 0., 10. );
    _hEdepROMC = tfs->make<TH1F>( "hEdepROMC", "Direct energy deposition in the APD", 100, 0., 10. );


    _hRCEdep    = tfs->make<TH1F>( "hRCEdep",
                                   "Total energy deposition in calorimeter reconstructed from ROs",
                                   2400, 0., 2400. );

    _hRCTime    = tfs->make<TH1F>( "hRCTime",
                                   "Hit time in calorimeter reconstructed from ROs",
                                   2400, 0., 2400. );

    _hRCNCrystals = tfs->make<TH1F>( "hRCNCrystals",
                                     "Number of crystals reconstructed from ROs",
                                     50, 0., 50. );

    _hRCEdepMC  = tfs->make<TH1F>( "hRCEdepMC",
                                   "Total energy deposition in calorimeter reconstructed from raw ROs MC",
                                   240, 0., 240. );

    _hRCTimeMC  = tfs->make<TH1F>( "hRCTimeMC",
                                   "Hit time in calorimeter reconstructed from raw ROs MC",
                                   2400, 0., 2400. );

    _hRCNCrystalsMC = tfs->make<TH1F>( "hRCNCrystalsMC",
                                       "Number of crystals reconstructed from raw ROs MC",
                                       50, 0., 50. );

    // Stopping target histograms

    _hTargetEdep = tfs->make<TH1F>( "hTargetEdep",
                                    "Energy deposition in the stopping target",
                                    100, 0., 5. );
    _hTargetPathLength = tfs->make<TH1F>( "hTargetPathLength",
                                          "Path length in the stopping target",
                                          100, 0., 5. );
    _hTargetNfoils = tfs->make<TH1F>( "hTargetNfoils",
                                      "Number of stopping target foils crossed by particle",
                                      20, 0., 20. );
    _hTargetNfoils2D = tfs->make<TH2F>( "hTargetNfoils2D",
                                        "Number of stopping target foils vs foil of origin",
                                        20, 0., 20., 20, 0, 20. );

    // Create tracker ntuple.
    _ntup = tfs->make<TNtuple>( "ntup", "Hit ntuple",
                                          "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec:lay:pdgId:genId:edep:p:step:hwz");

    // Create a TGraph;
    // - Syntax to set name and title is weird; that's just root.
    // - Must append to the output file by hand.
    _xyHits = tfs->make<TGraph>(_xyHitsMax);
    _xyHits->SetName("xyHits");
    _xyHits->SetTitle("Y vs X for StepPointMC");
    gDirectory->Append(_xyHits);

    // CRV

    _hCRVMultiplicity = tfs->make<TH1F>( "hCRVMultiplicity", "CRV StepPointMC per Bar", 5000,  0.,  5000. );

    // Create CRV ntuple.
    _ntupCRV = tfs->make<TNtuple>( "ntupCRV", "CRV Hit ntuple",
                                          "evt:trk:bid:hx:hy:hz:bx:by:bz:dx:dy:dz:lx:ly:lz:time:shld:mod:lay:pdgId:genId:edep:p:step");

  }

  void ReadBack::analyze(const art::Event& event) {

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Inquire about the completion status of G4.
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    if ( _nAnalyzed < _maxFullPrint ){
      cerr << g4Status << endl;
    }

    // Abort if G4 did not complete correctly.
    // Use your own judgement about whether to abort or to continue.
    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting ReadBack::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    // Call code appropriate for the tracker that is installed in this job.
    art::ServiceHandle<GeometryService> geom;
    if( geom->hasElement<LTracker>() || geom->hasElement<TTracker>() ){
      doLTracker(event);
    }
    else if ( geom->hasElement<ITracker>() ){
      doITracker(event);
    }

    doCalorimeter(event);

    doStoppingTarget(event);

    if( geom->hasElement<CosmicRayShield>() ) {
      doCRV(event);
    }
  }

  void ReadBack::doCalorimeter(const art::Event& event) {

    // Get handles to calorimeter collections
    art::Handle<CaloHitCollection> caloHits;
    art::Handle<CaloHitMCTruthCollection> caloMC;
    event.getByType(caloHits);
    event.getByType(caloMC);

    // Find pointers to the original G4 steps
    art::Handle<DPIndexVectorCollection> crystalPtr;
    art::Handle<DPIndexVectorCollection> readoutPtr;

    // One can use simple approach - find collection directly by label, e.g.:
    // event.getByLabel("CaloReadoutHitsMaker","CaloHitMCCrystalPtr",crystalPtr);
    // event.getByLabel("CaloReadoutHitsMaker","CaloHitMCReadoutPtr",readoutPtr);
    // The following code shows how to find collection only by name, not knowing
    // the producer module name

    vector<art::Handle<DPIndexVectorCollection> > ptr_coll;
    event.getManyByType(ptr_coll);
    for( size_t i=0; i<ptr_coll.size(); ++i ) {
      if(ptr_coll[i].provenance()->productInstanceName()=="CaloHitMCCrystalPtr") crystalPtr = ptr_coll[i];
      if(ptr_coll[i].provenance()->productInstanceName()=="CaloHitMCReadoutPtr") readoutPtr = ptr_coll[i];
    }

    // Find original G4 steps in the APDs
    art::Handle<StepPointMCCollection> rohits;
    event.getByLabel(_g4ModuleLabel,"calorimeterRO",rohits);

    bool haveCalo = ( caloHits.isValid() && caloMC.isValid() );

    if( ! haveCalo) return;

    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    double totalEdep = 0.0;
    double simEdep = 0.0;
    map<int,int> hit_crystals;

    for ( size_t i=0; i<caloHits->size(); ++i ) {
      totalEdep += caloHits->at(i).energyDep();
      simEdep += caloMC->at(i).energyDep();

      int roid = caloHits->at(i).id();
      int cid = cg->getCrystalByRO(roid);
      hit_crystals[cid] = 1;
    }

    _hEdep->Fill(totalEdep);
    _hEdepMC->Fill(simEdep);
    _hNcrystal->Fill(hit_crystals.size());

    // Fill number of G4 steps in crystal and APD, per APD
    if( crystalPtr.isValid() && readoutPtr.isValid() ) {
      for ( size_t i=0; i<caloHits->size(); ++i ) {
        _hNcrstep->Fill(crystalPtr->at(i).size());
        _hNrostep->Fill(readoutPtr->at(i).size());
      }
    }

    // Calculate direct energy deposition in the APD
    // This is an example how one can propagate back to the original
    // G4 data
    if( readoutPtr.isValid() && rohits.isValid() ) {
      for ( size_t i=0; i<caloHits->size(); ++i ) {
        // Get vector of pointer to G4 steps in APDs for calorimeter hit #i
        const DPIndexVector & ptr = readoutPtr->at(i);
        // Skip calorimeter hits without G4 step in APD (for these hit
        // no charged particle crossed APD)
        if( ptr.size()<=0 ) continue;
        // Accumulator to count total energy deposition
        double esum = 0;
        // Loop over that vector, get each G4 step and accumulate energy deposition
        for( size_t j=0; j<ptr.size(); j++ ) {
          const StepPointMC & rohit = rohits->at(ptr[j].index);
          esum += rohit.eDep();
        }
        // Fill histogram
        _hEdepROMC->Fill(esum);
      }
    }

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<caloHits->size(); ++i ) {
        CaloHit const & hit = (*caloHits).at(i);
        cout << "Readback: " << hit << endl;
      }
    }

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<caloMC->size(); ++i ) {
        CaloHitMCTruth const & hit = (*caloMC).at(i);
        cout << "Readback: " << hit << endl;
      }
    }


    // caloCrystalHits
    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;

    event.getByType(caloCrystalHits);
    if (!caloCrystalHits.isValid()) {
      _diagLevel > 0 && cout << __func__ << ": NO CaloCrystalHits" << endl;
      return;
    }

    totalEdep = 0.0;
    simEdep = 0.0;

    typedef multimap<int,size_t> hitCrystalsMultiMap;
    hitCrystalsMultiMap hitCrystals;

    _diagLevel > 0 &&
      cout << __func__ << ": caloCrystalHits->size() " << caloCrystalHits->size() << endl;

    for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {

      CaloCrystalHit const & hit = (*caloCrystalHits).at(i);

      totalEdep += hit.energyDep();
      _diagLevel > 0 && cout << __func__ << ": (*caloCrystalHits)[i].id(): "
                             << hit.id() << endl;

      // check if the crystal is there already (it may be ok if the timing is different)

      hitCrystalsMultiMap::const_iterator pos = hitCrystals.find(hit.id());

      if ( pos != hitCrystals.end() ) {

        _diagLevel > 0 && cout << __func__ << ": Already saw "
                               << (*caloCrystalHits).at(pos->second) << endl;

      }

      _diagLevel > 0 && cout << __func__ << ": Inserting   " << hit << endl;

      hitCrystals.insert(pair<int,size_t>(hit.id(),i));
      _hRCTime->Fill(hit.time());

    }

    _diagLevel > 0 && cout << __func__ << ": hitCrystals.size()     "
                           << hitCrystals.size() << endl;

    _hRCEdep->Fill(totalEdep);
    _hRCNCrystals->Fill(hitCrystals.size());

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<caloCrystalHits->size(); ++i ) {
        CaloCrystalHit const & hit = (*caloCrystalHits).at(i);
        cout << "Readback: " << hit << endl;
      }
    }

    // templetize it?
    // caloCrystalOnlyHits
    art::Handle<CaloCrystalOnlyHitCollection>  caloCrystalOnlyHits;

    event.getByType(caloCrystalOnlyHits);
    if (!caloCrystalOnlyHits.isValid()) {
      _diagLevel > 0 && cout << __func__ << ": NO CaloCrystalOnlyHits" << endl;
      return;
    }

    simEdep = 0.0;

    hitCrystals.clear();

    _diagLevel > 0 &&
      cout << __func__ << ": caloCrystalOnlyHits->size() " << caloCrystalOnlyHits->size() << endl;

    for ( size_t i=0; i<caloCrystalOnlyHits->size(); ++i ) {

      CaloCrystalOnlyHit const & hit = (*caloCrystalOnlyHits).at(i);

      simEdep += hit.energyDep();
      _diagLevel > 0 && cout << __func__ << ": (*caloCrystalOnlyHits)[i].id(): "
                             << hit.id() << endl;

      // check if the crystal is there already (it may be ok if the timing is different)

      hitCrystalsMultiMap::const_iterator pos = hitCrystals.find(hit.id());

      if ( pos != hitCrystals.end() ) {

        _diagLevel > 0 && cout << __func__ << ": Already saw "
                               << (*caloCrystalOnlyHits).at(pos->second) << endl;

      }

      _diagLevel > 0 && cout << __func__ << ": Inserting   " << hit << endl;

      hitCrystals.insert(pair<int,size_t>(hit.id(),i));
      _hRCTimeMC->Fill(hit.time());

    }

    _diagLevel > 0 && cout << __func__ << ": hitCrystals.size()     "
                           << hitCrystals.size() << endl;

    _hRCEdepMC->Fill(simEdep);
    _hRCNCrystalsMC->Fill(hitCrystals.size());

    if ( _diagLevel > -1 && _nAnalyzed < _maxFullPrint ){
      for ( size_t i=0; i<caloCrystalOnlyHits->size(); ++i ) {
        CaloCrystalOnlyHit const & hit = (*caloCrystalOnlyHits).at(i);
        cout << "Readback: " << hit << endl;
      }
    }

  }

  void ReadBack::doLTracker(const art::Event& event){

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    // Get handles to the generated and simulated particles.
    art::Handle<ToyGenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel,genParticles);

    art::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // A silly example just to show that we have a messsage logger.
    if ( hits->size() > 300 ){
      mf::LogWarning("HitInfo")
        << "Number of hits "
        << hits->size()
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cet::exception("RANGE")
        << "Way too many hits in this event.  Something is really wrong."
        << hits->size();
    }

    // ntuple buffer.
    float nt[_ntup->GetNvar()];

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get the straw information:
      const Straw&      straw = tracker.getStraw( hit.strawIndex() );
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();

      // Count how many nearest neighbours are also hit.
      int nNeighbours = countHitNeighbours( straw, hits );

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the straw.  Should be 2.5 mm.
      double s = w.dot(pos-mid);
      CLHEP::Hep3Vector point = pos - (mid + s*w);

      // The simulated particle that made this hit.
      SimParticleCollection::key_type trackId(hit.trackId());

      // Debug info

      StrawDetail const& strawDetail = straw.getDetail();

//       // remove this for production, intended for transportOnly.py
//       if (pca.dca()>strawDetail.innerRadius() || abs(point.mag()- strawDetail.innerRadius())>1.e-6 ) {

//         cerr << "*** Bad hit?: "
//              << event.id().event() << " "
//              << i                  <<  " "
//              << hit.trackId()      << "   "
//              << hit.volumeId()     << " "
//              << straw.Id()         << " | "
//              << pca.dca()          << " "
//              << pos                << " "
//              << mom                << " "
//              << point.mag()        << " "
//              << hit.eDep()         << " "
//              << s
//              << endl;
//       }

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
          genId = gen.generatorId();
        }
      }

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hEnergyDep->Fill(hit.eDep()/keV);
      _hTime->Fill(hit.time());
      _hHitNeighbours->Fill(nNeighbours);
      _hCheckPointRadius->Fill(point.mag());
      _hCheckPointRadiusW->Fill(point.mag()/strawDetail.innerRadius());
      _hCheckPointWireZ->Fill(s/straw.getHalfLength());

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());
      _hDriftDistW->Fill(pca.dca()/strawDetail.innerRadius());

      _hStepLength->Fill( hit.stepLength() );

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = hit.trackId().asInt();
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      nt[6]  = mid.x();
      nt[7]  = mid.y();
      nt[8]  = mid.z();
      nt[9]  = pca.dca();
      nt[10] = hit.time();
      nt[11] = straw.Id().getDevice();
      nt[12] = straw.Id().getSector();
      nt[13] = straw.Id().getLayer();
      nt[14] = pdgId;
      nt[15] = genId.Id();
      nt[16] = hit.eDep()/keV;
      nt[17] = mom.mag();
      nt[18] = hit.stepLength();
      nt[19] = s/straw.getHalfLength();

      _ntup->Fill(nt);

      // Fill the TGraph; need to manage size by hand.
      if ( _xyHitCount < _xyHitsMax ){
        _xyHits->SetPoint( _xyHitCount++, pos.x(), pos.y() );
      }

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){

        cerr << "Readback"
             << " hit: "
             << event.id().event() << " "
             << i                  <<  " "
             << hit.trackId()      << "   "
             << hit.volumeId()     << " "
             << straw.Id()         << " | "
             << pca.dca()          << " "
             << pos                << " "
             << mom                << " "
             << point.mag()        << " "
             << hit.eDep()         << " "
             << s
             << endl;
      }

    } // end loop over hits.


    // Additional printout and histograms about the simulated particles.
    if ( haveSimPart && (_nAnalyzed < _maxFullPrint) ){

      ConditionsHandle<ParticleDataTable> pdt("ignored");

      for ( SimParticleCollection::const_iterator i=simParticles->begin();
            i!=simParticles->end(); ++i ){

        SimParticle const& sim = i->second;

        if ( sim.madeInG4() ) {

          _hMomentumG4->Fill( sim.startMomentum().rho() );

        } else {

          // Particle Data Group Id number of this SimParticle
          int pdgId = sim.pdgId();

          // Name of this particle type.
          ParticleDataTable::maybe_ref particle = pdt->particle(pdgId);
          string pname = particle ?
            particle.ref().name() :
            unknownPDGIdName(pdgId);

          // Information about generated particle.
          ToyGenParticle const& gen = genParticles->at(sim.generatorIndex());
          GenId genId(gen.generatorId());

          // Physical volume in which this track started.
          PhysicalVolumeInfo const& volInfo = volumes->at(sim.startVolumeIndex());
          PhysicalVolumeInfo const& endInfo = volumes->at(sim.endVolumeIndex());

          cerr << "Readback"
               << " Simulated Particle: "
               << i->first            << " "
               << pdgId               << " "
               << genId.name()        << " "
               << sim.startPosition() << " "
               << volInfo.name()      << " "
               << volInfo.copyNo()    << " | "
               << endInfo.name()      << " "
               << endInfo.copyNo()    << " | "
               << sim.stoppingCode()
               << endl;
        }

      }
    }

  } // end doLTracker

  void ReadBack::doITracker(const art::Event& event){

    // Instance name of the module that created the hits of interest;
    //static const string creatorName("g4run");

    // Gometry for the ITracker.
    GeomHandle<ITracker> itracker;
    CellGeometryHandle *itwp = itracker->getCellGeometryHandle();

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    // Get handles to the generated and simulated particles.
    art::Handle<ToyGenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel,genParticles);

    art::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // A silly example just to show that we have a messsage logger.
    if ( hits->size() > 300 ){
      mf::LogWarning("HitInfo")
        << "Number of hits "
        << hits->size()
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cet::exception("RANGE")
        << "Way too many hits in this event.  Something is really wrong."
        << hits->size();
    }

    // ntuple buffer.
    //float nt[13];
    float nt[_ntup->GetNvar()];

    // Loop over all hits.
    int n(0);
    StepPointMCCollection::const_iterator i = hits->begin();
    StepPointMCCollection::const_iterator e = hits->end();
    for ( ; i!=e; ++i){

      // Aliases, used for readability.
      const StepPointMC& hit = *i;

      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;

      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      itwp->SelectCellDet(hit.volumeId());
      //Get the cell information.
          boost::shared_ptr<mu2e::Cell> cell = itwp->GetITCell();
          //Cell const& cell = itwp->GetITCell();
          CLHEP::Hep3Vector mid = cell->getMidPoint();
          CLHEP::Hep3Vector w   = cell->getDirection();

      // Count how many nearest neighbours are also hit.
      //    int nNeighbours = countHitNeighbours( cell, hits );

      // Compute an estimate of the drift distance.
         TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the cell.  Should be 2.5 mm.
          double s = w.dot(pos-mid);
          CLHEP::Hep3Vector point = pos - (mid + s*w);

          // The simulated particle that made this hit.
          SimParticleCollection::key_type trackId(hit.trackId());

      // I don't understand the distribution of the time variable.
      // I want it to be the time from the start of the spill.
      // It appears to be the time since start of tracking.

      // Fill some histograms
      //    _hRadius->Fill(pos.perp());
      _hTime->Fill(hit.time());
      //    _hHitNeighbours->Fill(nNeighbours);
      //    _hCheckPointRadius->Fill(point.mag());

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      //    _hDriftDist->Fill(pca.dca());
      double distUnit = (itracker->isExternal()) ? 1.0*CLHEP::cm : 1.0*CLHEP::mm ;
      double invDistUnit = 1.0/distUnit;
      double hitpos[3] = {pos.x()*invDistUnit,pos.y()*invDistUnit,pos.z()*invDistUnit};
      double lclhitpos[3] = {0.0,0.0,0.0};
      itwp->Global2Local(hitpos,lclhitpos);

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
          genId = gen.generatorId();
        }
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = hit.trackId().asInt();
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      nt[6]  = distUnit*lclhitpos[0];//mid.x();
      nt[7]  = distUnit*lclhitpos[1];//mid.y();
      nt[8]  = distUnit*lclhitpos[2];//mid.z();
      nt[9]  = distUnit*itwp->DistFromWire(hitpos);//pca.dca();
      nt[10] = hit.time();
      nt[11] = itwp->GetITCell()->Id().getCell();
      nt[12] = itwp->GetITCell()->Id().getLayer();
      nt[13] = itwp->GetITCell()->Id().getLayerId().getSuperLayer();
      nt[14] = pdgId;
      nt[15] = genId.Id();
      nt[16] = hit.eDep()/keV;
      nt[17] = mom.mag();
      nt[18] = hit.stepLength();
      nt[19] = s/itwp->GetITCell()->getHalfLength();

      //    if (nt[6]<-3.0) {
      //
      //    }

      _ntup->Fill(nt);

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){
        cerr << "ReadBack"
             << " hit: "
             << event.id().event() << " "
             << n++ <<  " "
             << hit.trackId()  << "   "
             << hit.volumeId() << " | "
          /*<< pca.dca()   << " "*/
             << pos  << " "
             << mom  << " "
          /*<< point.mag()*/
             << endl;
      }

    } // end loop over hits.

  }  // end doITracker


  // Count how many of this straw's nearest neighbours are hit.
  // If we have enough hits per event, it will make sense to make
  // a class to let us direct index into a list of which straws have hits.
  int ReadBack::countHitNeighbours( Straw const& straw,
                                    art::Handle<StepPointMCCollection>& hits ){

    int count(0);
    vector<StrawIndex> const& nearest = straw.nearestNeighboursByIndex();
    for ( vector<int>::size_type ihit = 0;
          ihit<nearest.size(); ++ihit ){

      StrawIndex idx = nearest[ihit];

      for( StepPointMCCollection::const_iterator
             i = hits->begin(),
             e = hits->end(); i!=e ; ++i ) {
        const StepPointMC& hit = *i;
        if ( hit.strawIndex() == idx ){
          ++count;
          break;
        }
      }

    }
    return count;
  }  // end countHitNeighbours

  //
  // Example of how to read information about stopping target
  //
  // Here we assume that the primary particles are the conversion
  // electrons generated in the stopping target.
  //
  void ReadBack::doStoppingTarget(const art::Event& event) {

    // Find original G4 steps in the stopping target
    art::Handle<StepPointMCCollection> sthits;
    event.getByLabel(_g4ModuleLabel,_targetStepPoints,sthits);

    // SimParticles container
    art::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);
    if( !(simParticles.isValid()) || simParticles->empty() ) return;

    // Loop over all hits in the stopping target. Check, that the
    // hit belongs to primary particle. If so, calculate total energy
    // deposition in the target, total path length and the number of
    // foils, crossed by the electron.

    int id_start = -1;
    double eDep=0.0, pathLength=0.0;
    map<int,int> foils;

    // Loop over all hits in the stopping target
    for ( size_t i=0; i<sthits->size(); ++i ){

      // This is G4 hit (step) in the target
      const StepPointMC& hit = (*sthits)[i];

      // Here we select only those hits, which are generated by
      // the primary track - which is assumed to be the original
      // particle, generated by ConversionGun
      SimParticleCollection::key_type trackId = hit.trackId();
      if( trackId.asInt() != 1 ) continue;

      // Here we require that there is information about the primary
      // particle in the SimParticle collection. It is not neccessary for
      // this example, but it is typical requirement in the real analysis
      SimParticle const* sim = simParticles->findOrNull(trackId);
      if( !sim ) continue;

      // Get the foil id where the hit occured. If it is the first hit,
      // remember this id as the source foil id.
      int id = hit.volumeId();
      if( id_start<0 ) id_start=id;

      // Here we calculate number of steps in each foil. There could be
      // many hits in each foil. This number is not used in the example,
      // it is calculated here just as an example. But we use this map to
      // calculate number of foils with the hits. Be aware, that if particle
      // crosses foil without energy deposition (just passes it), it still
      // counts as a hit. Here we record all foils particle crosses.
      // If we want to record only those foils where particle had
      // interactions, we would need to do the following:
      //      if( hit.totalEDep()>0 ) foils[id]++;
      foils[id]++;

      // Calculate total energy loss and path length primary particle
      // has in the target.
      eDep += hit.totalEDep();
      pathLength += hit.stepLength();

    }

    // Number of crossed foils
    int nfoil = foils.size();

    // Fill histograms
    if( id_start>=0 ) {
      _hTargetEdep->Fill(eDep);
      _hTargetPathLength->Fill(pathLength);
      _hTargetNfoils->Fill(nfoil);
      _hTargetNfoils2D->Fill(id_start,nfoil);
    }

  } // end doStoppingTarget

  void ReadBack::endJob(){
    cout << "ReadBack::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
         << endl;
  }

  void ReadBack::doCRV(const art::Event& event){

    // Get a reference to CosmicRayShield (it contains crv)

    GeomHandle<CosmicRayShield> cosmicRayShieldGeomHandle;

    std::vector<CRSScintillatorBar> const & allBars =
      cosmicRayShieldGeomHandle->getAllCRSScintillatorBars();

    CRSScintillatorBarDetail const & barDetail =
      cosmicRayShieldGeomHandle->getCRSScintillatorBarDetail();

    CLHEP::Hep3Vector barLengths = CLHEP::Hep3Vector(barDetail.getHalfLengths()[0],
                                                     barDetail.getHalfLengths()[1],
                                                     barDetail.getHalfLengths()[2]);

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_crvStepPoints,hits);

    // Get handles to the generated and simulated particles.
    art::Handle<ToyGenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel,genParticles);

    art::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    // A silly example just to show that we have a messsage logger.
    if ( hits->size() > 300 ){
      mf::LogWarning("HitInfo")
        << "Number of CRV hits "
        << hits->size()
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cet::exception("RANGE")
        << "Way too many CRV hits in this event.  Something is really wrong."
        << hits->size();
    }

    if ( _nAnalyzed < _maxFullPrint ){
      cerr << g4Status << endl;
    }

    // ntuple buffer.
    float nt[_ntupCRV->GetNvar()];

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get the CRSScintillatorBar information:
      const CRSScintillatorBar&  bar = allBars.at( hit.volumeId() );
      CLHEP::Hep3Vector const &  mid = bar.getGlobalOffset();

      CLHEP::HepRotationX RX(bar.getGlobalRotationAngles()[0]);
      CLHEP::HepRotationY RY(bar.getGlobalRotationAngles()[1]);
      CLHEP::HepRotationZ RZ(bar.getGlobalRotationAngles()[2]);

      CLHEP::HepRotation barInvRotation((RX*RY*RZ).inverse());
      CLHEP::Hep3Vector hitLocal  = barInvRotation*(pos-mid);
      CLHEP::Hep3Vector hitLocalN = CLHEP::Hep3Vector(hitLocal.x()/barLengths.x(),
                                                      hitLocal.y()/barLengths.y(),
                                                      hitLocal.z()/barLengths.z());

      // The simulated particle that made this hit.
      SimParticleCollection::key_type trackId(hit.trackId());

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
            genId = gen.generatorId();
          }

      }

      // Fill the ntuple.
      nt[ 0] = event.id().event();
      nt[ 1] = hit.trackId().asInt();
      nt[ 2] = hit.volumeId();
      nt[ 3] = pos.x();
      nt[ 4] = pos.y();
      nt[ 5] = pos.z();
      nt[ 6] = mid.x();
      nt[ 7] = mid.y();
      nt[ 8] = mid.z();
      nt[ 9] = fabs((pos-mid).x());
      nt[10] = fabs((pos-mid).y());
      nt[11] = fabs((pos-mid).z());
      nt[12] = hitLocalN.x();
      nt[13] = hitLocalN.y();
      nt[14] = hitLocalN.z();
      nt[15] = hit.time();
      nt[16] = bar.Id().getShieldNumber();
      nt[17] = bar.Id().getModuleNumber();
      nt[18] = bar.Id().getLayerNumber();
      nt[19] = pdgId;
      nt[20] = genId.Id();
      nt[21] = hit.eDep()/keV;
      nt[22] = mom.mag();
      nt[23] = hit.stepLength();

      _ntupCRV->Fill(nt);

      _hCRVMultiplicity->Fill(hit.volumeId());

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){

        cerr << "Readback"
             << " hit: "
             << event.id().event() << " "
             << i                  <<  " "
             << hit.trackId()      << "   "
             << hit.volumeId()     << " "
             << bar.Id()           << " | "
             << pos                << " "
             << mid                << " "
             << (mid-pos)          << " | "
             << hitLocalN          << " | "
             << mom                << " "
             << hit.totalEDep()/keV << " "
             << hit.stepLength()
             << endl;
      }

    } // end loop over hits.

  } // end doCRV

}  // end namespace mu2e
