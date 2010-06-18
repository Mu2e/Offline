//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: ReadBack.cc,v 1.9 2010/06/18 18:45:57 genser Exp $
// $Author: genser $
// $Date: 2010/06/18 18:45:57 $
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
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"

// Root includes.
#include "TH1F.h"
#include "TNtuple.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  ReadBack::ReadBack(edm::ParameterSet const& pset) : 
    _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
    _minimumEnergy(pset.getParameter<double>("minimumEnergy")),
    _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
    _nAnalyzed(0),
    _hRadius(0),
    _hTime(0),
    _hMultiplicity(0),
    _hDriftDist(0),
    _hxHit(0),
    _hyHit(0),
    _hzHit(0),
    _hHitNeighbours(0),
    _hCheckPointRadius(0),
    _ntup(0){
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

    _hxHit         = tfs->make<TH1F>( "hxHit",  "X of Hit;(mm)",                  100, -1000., 1000. );
    _hyHit         = tfs->make<TH1F>( "hyHit",  "Y of Hit;(mm)",                  100, -1000., 1000. );
    _hzHit         = tfs->make<TH1F>( "hzHit",  "Z of Hit;(mm)",                  100, -1400., 1400. );

    _hHitNeighbours    = tfs->make<TH1F>( "hHitNeighbours",  "Number of hit neighbours",
                                          10, 0., 10. );

    _hCheckPointRadius = tfs->make<TH1F>( "hCheckPointRadius",  "Radius of Reference point; (mm)",
                                          100, 2.25, 2.75 );

    _hMomentumG4 = tfs->make<TH1F>( "hMomentumG4",  "Mommenta of particles created inside G4; (MeV)",
                                    100, 0., 100. );
    _hStepLength = tfs->make<TH1F>( "hStepLength",  "G4 Step Length in Sensitive Detector; (mm)",
                                    100, 0., 10. );

    // Create an ntuple.
    _ntup           = tfs->make<TNtuple>( "ntup", "Hit ntuple", 
                      "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec:pdgId:genId:edep:p:step");
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

  }

  void ReadBack::doLTracker(const edm::Event& event){

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,hits);

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
    if ( hits->size() > 75 ){
      edm::LogWarning("HitInfo")
        << "Number of hits "
        << hits->size() 
        << " is too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cms::Exception("RANGE")
        << "Way too many hits in this event.  Something is really wrong."
        << hits->size();
    }

    // ntuple buffer.
    float nt[18];

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
      int trackId = hit.trackId();

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
      
      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());

      _hStepLength->Fill( hit.stepLength() );

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = hit.trackId();
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
      nt[13] = pdgId;
      nt[14] = genId.Id();
      nt[15] = hit.eDep()/keV;
      nt[16] = mom.mag();
      nt[17] = hit.stepLength();

      _ntup->Fill(nt);

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){

        cerr << "Readback hit: "
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

      for ( int i=0; i<simParticles->size(); ++ i){

        SimParticle const& sim = simParticles->at(i);

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

          cerr << "Simulated Particle: " 
               << i                   << " "
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
    event.getByLabel(creatorName,hits);

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
      nt[1]  = hit.trackId();
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
        cerr << "Readback hit: "
             << event.id().event() << " "
             << n++ <<  " "
             << hit.trackId()    << "   "
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
