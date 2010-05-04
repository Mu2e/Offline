//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CosmicTuple.cc,v 1.1 2010/05/04 01:33:46 yury Exp $
// $Author: yury $
// $Date: 2010/05/04 01:33:46 $
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
#include "Mu2eG4/src/CosmicTuple.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

// Root includes.
#include "TH1F.h"
#include "TNtuple.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  CosmicTuple::CosmicTuple(edm::ParameterSet const& pset) : 
    _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
    _minimumEnergy(pset.getParameter<double>("minimumEnergy")),
    _nAnalyzed(0),
    _ntupTrk(0)
  {
  }
  
  void CosmicTuple::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;
    
    // Create an ntuple.
    _ntupTrk = tfs->make<TNtuple>( "ntupTrk", "Trk ntuple", 
				   "evt:trk:pid:px:py:pz:pmag:genId:pidGen:eGen:thGen:xGen:yGen:zGen:nHits:hxMin:hyMin:hzMin:hxMax:hyMax:hzMax:tMin:tMax");

  }

  void CosmicTuple::analyze(const edm::Event& event, edm::EventSetup const&) {
    
    // Maintain a counter for number of events seen.
    ++_nAnalyzed;
    
    if ( _nAnalyzed % 10000 == 0 ) {
      edm::LogInfo("CosmicTuple")
        << "Processing event " << _nAnalyzed;
    }

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,hits);

    // Get handles to the generated and simulated particles.
    edm::Handle<ToyGenParticleCollection> genParticles;
    event.getByType(genParticles);

    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() );
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cms::Exception("RANGE")
        << "Way too many hits in this event.  Something is really wrong."
        << hits->size();
    }


    // ntuple buffer.
    float ntT[23];   

    int oldTrk = -1;
    bool first = true;
    Hep3Vector ptrk(0,0,0);
    int nHits = 0;
    float hxMin = 1e6;
    float hyMin = 1e6;
    float hzMin = 1e6;
    float tMin = 1e6;
    float hxMax = -1e6;
    float hyMax = -1e6;
    float hzMax = -1e6;
    float tMax = -1e6;
 
    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){
      
      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;
      
      // Get the hit information.
      const Hep3Vector& pos = hit.position();
      const Hep3Vector& mom = hit.momentum();
      
      // The simulated particle that made this hit.
      int trackId = hit.trackId();

      if ( oldTrk != trackId || i == hits->size()-1 ) {
	// special case
	if ( i == hits->size()-1 ) {
	  oldTrk = trackId;
	  if ( pos.x() < hxMin ) hxMin = pos.x();
	  if ( pos.y() < hyMin ) hyMin = pos.y();
	  if ( pos.z() < hzMin ) hzMin = pos.z();
	  if ( hit.time() < tMin ) tMin = hit.time();
	  if ( pos.x() > hxMax ) hxMax = pos.x();
	  if ( pos.y() > hyMax ) hyMax = pos.y();
	  if ( pos.z() > hzMax ) hzMax = pos.z();
	  if ( hit.time() > tMax ) tMax = hit.time();
	  nHits++;
	}

        if ( !first ) {
          // dump previous track

	  // Default values for these, in case information is not available.
	  int pdgId = -1;
	  float eGen = 0;
	  float thGen = 0;
	  Hep3Vector posGen(0,0,0);
	  GenId idGen;
	  int pidGen = -1;

	  if ( haveSimPart && oldTrk >= 0 && oldTrk < simParticles->size() ){
	    const SimParticle* sim = &(simParticles->at(oldTrk));

	    // PDG Particle Id of the sim particle that made this hit.
	    pdgId = sim->pdgId();
      
	    // find the generated parent
	    const SimParticle* p = sim;
	    int gTrk = p->generatorIndex();
	    int depth = 0;
	    while ( gTrk < 0 && depth < 100 ) {
	      if ( p->hasParent() ) {
		int pId = p->parentId();
		if ( pId < 0 || pId >= simParticles->size() ) break;
		p = &(simParticles->at(pId));
		gTrk = p->generatorIndex();
		depth++;
	      } else {
		break;
	      }
	    }

	    // store generator info
	    if ( gTrk >= 0 && gTrk < genParticles->size() ) {
	      ToyGenParticle const& genpart = genParticles->at(gTrk);
	      idGen = genpart._generatorId;
	      pidGen = genpart._pdgId;
	      HepLorentzVector p4gen = genpart._momentum;
	      Hep3Vector y(0,-1,0);
	      eGen = p4gen.e();
	      thGen = y.angle(p4gen.vect());
	      posGen = genpart._position;
	    }

	  }

	  ntT[0]  = event.id().event();
          ntT[1]  = oldTrk;
	  ntT[2]  = pdgId;
          ntT[3]  = ptrk.x();
          ntT[4]  = ptrk.y();
          ntT[5]  = ptrk.z();
          ntT[6]  = ptrk.mag();
          ntT[7]  = idGen.Id();
          ntT[8]  = pidGen;
          ntT[9]  = eGen;
          ntT[10] = thGen;
          ntT[11] = posGen.x();
          ntT[12] = posGen.y();
          ntT[13] = posGen.z();
          ntT[14] = nHits;
          ntT[15] = hxMin;
          ntT[16] = hyMin;
          ntT[17] = hzMin;
          ntT[18] = hxMax;
          ntT[19] = hyMax;
          ntT[20] = hzMax;
          ntT[21] = tMin;
          ntT[22] = tMax;

          _ntupTrk->Fill(ntT);
        }

        oldTrk = hit.trackId();
        ptrk = mom;
        nHits = 0;
        hxMin = 1e6;
        hyMin = 1e6;
        hzMin = 1e6;
        tMin = 1e6;
        hxMax = -1e6;
        hyMax = -1e6;
        hzMax = -1e6;
        tMax=  -1e6;
      }

      if ( pos.x() < hxMin ) hxMin = pos.x();
      if ( pos.y() < hyMin ) hyMin = pos.y();
      if ( pos.z() < hzMin ) hzMin = pos.z();
      if ( hit.time() < tMin ) tMin = hit.time();
      if ( pos.x() > hxMax ) hxMax = pos.x();
      if ( pos.y() > hyMax ) hyMax = pos.y();
      if ( pos.z() > hzMax ) hzMax = pos.z();
      if ( hit.time() > tMax ) tMax = hit.time();
      nHits++;
      first=false;

    } // end loop over hits.


  } // end analyze

  
}  // end namespace mu2e
