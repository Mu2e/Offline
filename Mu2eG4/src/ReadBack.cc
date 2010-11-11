//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: ReadBack.cc,v 1.21 2010/11/11 21:23:47 genser Exp $
// $Author: genser $
// $Date: 2010/11/11 21:23:47 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

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
#include "ToyDP/inc/CaloCrystalHitMCTruthCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

// Root includes.
#include "TDirectory.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TGraph.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  ReadBack::ReadBack(edm::ParameterSet const& pset) : 

    // Run time parameters
    _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
    _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
    _caloCrystalHitsMaker(pset.getUntrackedParameter<string>("caloCrystalHitsMaker","CaloCrystalHitsMaker")),
    _minimumEnergy(pset.getParameter<double>("minimumEnergy")),
    _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
    _xyHitsMax(pset.getUntrackedParameter<int>("xyHitsMax",10000)),

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
    _hRCEdep(0),
    _hRCTime(0),
    _hRCNCrystals(0),
    _hRCEdepMC(0),
    _hRCTimeMC(0),
    _hRCNCrystalsMC(0),
    _ntup(0),
    _xyHits(0),

    // Remaining member data
    _xyHitCount(0){
  }
  
  void ReadBack::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;
    
    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius",       "Radius of Hits;(mm)",     100,  0., 1000. );
    _hEnergyDep    = tfs->make<TH1F>( "hEnergyDep",    "Energy Deposited;(keV)",  100,  0.,   10. );
    _hTime         = tfs->make<TH1F>( "hTime",         "Pulse Height;(ns)",       100,  0., 2000. );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event",          100,  0.,  100. );
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

    // Create an ntuple.
    _ntup           = tfs->make<TNtuple>( "ntup", "Hit ntuple", 
                      "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec:lay:pdgId:genId:edep:p:step:hwz");

    // Create a TGraph; 
    // - Syntax to set name and title is weird; that's just root.
    // - Must append to the output file by hand.
    _xyHits = tfs->make<TGraph>(_xyHitsMax);
    _xyHits->SetName("xyHits");
    _xyHits->SetTitle("Y vs X for StepPointMC");
    gDirectory->Append(_xyHits);

  }

  void ReadBack::analyze(const edm::Event& event, edm::EventSetup const&) {
    
    // Maintain a counter for number of events seen.
    ++_nAnalyzed;
    
    // Call code appropriate for the tracker that is installed in this job.
    edm::Service<GeometryService> geom;
    if( geom->hasElement<LTracker>() || geom->hasElement<TTracker>() ){
      doLTracker(event);
    }
    else if ( geom->hasElement<ITracker>() ){
      doITracker(event);
    }

    doCalorimeter(event);

  }

  void ReadBack::doCalorimeter(const edm::Event& event) {
    
    // Get handles to calorimeter collections
    edm::Handle<CaloHitCollection> caloHits;
    edm::Handle<CaloHitMCTruthCollection> caloMC;
    event.getByType(caloHits);
    event.getByType(caloMC);
    bool haveCalo = ( caloHits.isValid() && caloMC.isValid() );

    if( ! haveCalo) return;

    edm::Service<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    double totalEdep = 0.0;
    double simEdep = 0.0;
    map<int,int> hit_crystals;

    for ( size_t i=0; i<caloHits->size(); ++i ) {
      totalEdep += caloHits->at(i).energyDep();
      simEdep += caloMC->at(i).energyDep();

      int roid = caloHits->at(i).roId();
      int cid = cg->getCrystalByRO(roid);
      hit_crystals[cid] = 1;
    }

    _hEdep->Fill(totalEdep);
    _hEdepMC->Fill(simEdep);
    _hNcrystal->Fill(hit_crystals.size());


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
    edm::Handle<CaloCrystalHitCollection>  caloCrystalHits;

    event.getByType(caloCrystalHits);
    if (!caloCrystalHits.isValid()) {
      _diagLevel > 0 && cout << __func__ << ": NO CaloCrystalHits" << endl;
      return;
    }

    totalEdep = 0.0;
    simEdep = 0.0;

    typedef multimap<int,size_t> hitCrystalsMultiMap;
    multimap<int,size_t> hitCrystals;

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
    // caloCrystalHitMCTruths
    edm::Handle<CaloCrystalHitMCTruthCollection>  caloCrystalHitMCTruths;

    event.getByType(caloCrystalHitMCTruths);
    if (!caloCrystalHitMCTruths.isValid()) {
      _diagLevel > 0 && cout << __func__ << ": NO CaloCrystalHitMCTruths" << endl;
      return;
    }

    simEdep = 0.0;

    typedef multimap<int,size_t> hitCrystalsMultiMap;
    hitCrystals.clear();

    _diagLevel > 0 && 
      cout << __func__ << ": caloCrystalHitMCTruths->size() " << caloCrystalHitMCTruths->size() << endl;

    for ( size_t i=0; i<caloCrystalHitMCTruths->size(); ++i ) {

      CaloCrystalHitMCTruth const & hit = (*caloCrystalHitMCTruths).at(i);

      simEdep += hit.energyDep();
      _diagLevel > 0 && cout << __func__ << ": (*caloCrystalHitMCTruths)[i].id(): " 
                             << hit.id() << endl;

      // check if the crystal is there already (it may be ok if the timing is different)

      hitCrystalsMultiMap::const_iterator pos = hitCrystals.find(hit.id());

      if ( pos != hitCrystals.end() ) {

        _diagLevel > 0 && cout << __func__ << ": Already saw " 
               << (*caloCrystalHitMCTruths).at(pos->second) << endl;
        
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
      for ( size_t i=0; i<caloCrystalHitMCTruths->size(); ++i ) {
        CaloCrystalHitMCTruth const & hit = (*caloCrystalHitMCTruths).at(i);
        cout << "Readback: " << hit << endl;
      }
    }

  }

  void ReadBack::doLTracker(const edm::Event& event){

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    // Get handles to the generated and simulated particles.
    edm::Handle<ToyGenParticleCollection> genParticles;
    event.getByType(genParticles);

    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    edm::Handle<PhysicalVolumeInfoCollection> volumes;
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
      edm::LogWarning("HitInfo")
        << "Number of hits "
        << hits->size() 
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cms::Exception("RANGE")
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

//       // remove this for production, intended for transporOnly.py
//       if (pca.dca()>strawDetail.outerRadius() || abs(point.mag()- strawDetail.outerRadius())>1.e-6 ) {

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
          genId = gen._generatorId;
        }
      }

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hEnergyDep->Fill(hit.eDep()/keV);
      _hTime->Fill(hit.time());
      _hHitNeighbours->Fill(nNeighbours);
      _hCheckPointRadius->Fill(point.mag());
      _hCheckPointRadiusW->Fill(point.mag()/strawDetail.outerRadius());
      _hCheckPointWireZ->Fill(s/straw.getHalfLength());
      
      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());
      _hDriftDistW->Fill(pca.dca()/strawDetail.outerRadius());

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

          // Particle Data group Id number.
          int pdgId = sim.pdgId();
          
          // Information about generated particle.
          ToyGenParticle const& gen = genParticles->at(sim.generatorIndex());
          GenId genId(gen._generatorId);

          // Physical volume in which this track started.
          PhysicalVolumeInfo const& volInfo = volumes->at(sim.startVolumeIndex());

          cerr << "Readback"
               << " Simulated Particle: " 
               << i->first            << " "
               << pdgId               << " "
               << genId.name()        << " "
               << sim.startPosition() << " "
               << volInfo.name()      << " "
               << volInfo.copyNo()
               << endl;
        }

      }
    }

  } // end doLTracker

  void ReadBack::doITracker(const edm::Event& event){
     
    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");
 
    // Gometry for the ITracker.
    GeomHandle<ITracker> itracker;
    CellGeometryHandle *itwp = itracker->getCellGeometryHandle();
    
    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(creatorName,_trackerStepPoints,hits);

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // ntuple buffer.
    float nt[13];
    
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

      // Get the cell information.
      //    Cell const& cell = itracker->getCell( hit.volumeId() );
      //    CLHEP::Hep3Vector mid = cell.getMidPoint();
      //    CLHEP::Hep3Vector w   = cell.getDirection();

      // Count how many nearest neighbours are also hit.
      //    int nNeighbours = countHitNeighbours( cell, hits );

      // Compute an estimate of the drift distance.
      //   TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the cell.  Should be 2.5 mm.
      //    double s = w.dot(pos-mid);
      //    CLHEP::Hep3Vector point = pos - (mid + s*w);

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
      itwp->SelectWireDet(hit.volumeId());
      double distUnit = (itracker->isExternal()) ? 1.0*CLHEP::cm : 1.0*CLHEP::mm ;
      double invDistUnit = 1.0/distUnit;
      double hitpos[3] = {pos.x()*invDistUnit,pos.y()*invDistUnit,pos.z()*invDistUnit};
      double lclhitpos[3] = {0.0,0.0,0.0};
      itwp->Global2Local(hitpos,lclhitpos);
      
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
      nt[11] = 0;//cell.Id().getDevice();
      nt[12] = 0;//cell.Id().getSector();

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
                                    edm::Handle<StepPointMCCollection>& hits ){
    
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
  }
  
}  // end namespace mu2e
