//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.cc,v 1.10 2010/04/02 18:17:00 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/04/02 18:17:00 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <cassert>
#include <stdexcept>

// Framework includes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// G4 Includes
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

// Mu2e includes
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"

// ROOT includes
#include "TH1D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

  PrimaryGeneratorAction::PrimaryGeneratorAction( const string& generatorModuleLabel ):
    _generatorModuleLabel(generatorModuleLabel),
    _randomUnitSphere(0){

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    _particleDefinition = particleTable->FindParticle("chargedgeantino");

    // Generator for uniform coverage of a restricted region on a unit sphere.
    _randomUnitSphere =  auto_ptr<RandomUnitSphere>(new RandomUnitSphere ( -0.7, 0.7 ));

    // Book histograms.
    edm::Service<edm::TFileService> tfs;
    _totalMultiplicity = tfs->make<TH1D>( "totalMultiplicity", "Total Multiplicity", 20, 0, 20);

  }


  PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  }


  void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
  { 

    // For debugging.
    //testTrack(anEvent);

    fromEvent(anEvent);

  }

  // Copy generated particles from the event into G4.
  //
  // At the moment none of the supported generators make multi-particle vertices.
  // That may change in the future.
  // 
  void PrimaryGeneratorAction::fromEvent(G4Event* event){ 

    // Get the offsets to map from generator world to G4 world.
    G4ThreeVector const& cosmicReferencePlane = _world->getCosmicReferencePoint();
    G4ThreeVector const& detectorOrigin       = _world->getMu2eDetectorOrigin();
    G4ThreeVector const& primaryProtonGunOrigin      = _world->getPrimaryProtonGunOrigin();
    G4RotationMatrix const& primaryProtonGunRotation = _world->getPrimaryProtonGunRotation();

    // Get generated particles from the event.
    edm::Handle<ToyGenParticleCollection> handle;
    _event->getByLabel(_generatorModuleLabel,handle);

    // Fill multiplicity histogram.
    _totalMultiplicity->Fill( handle->size() );
  
    // For each generated particle, add it to the event.
    for ( unsigned int i=0; i<handle->size(); ++i){

      ToyGenParticle const& genpart = (*handle)[i];

      // Transform from generator coordinate system G4 world coordinate system.
      G4ThreeVector      pos(genpart._position);
      G4ThreeVector momentum(genpart._momentum.v());

      if( genpart._generatorId == GenId::conversionGun ||
          genpart._generatorId == GenId::particleGun ||
          genpart._generatorId == GenId::dio1 ||
          genpart._generatorId == GenId::ejectedProtonGun ||
          genpart._generatorId == GenId::pionCapture ||
          genpart._generatorId == GenId::piEplusNuGun){
        pos += detectorOrigin;
      } else if ( genpart._generatorId == GenId::cosmicToy ||
                  genpart._generatorId == GenId::cosmicDYB || 
                  genpart._generatorId == GenId::cosmic ){
        pos += cosmicReferencePlane;
      } else if ( genpart._generatorId == GenId::primaryProtonGun ){	
	//
	// need inverse; we define the physical object to point along its own z; 
	// then rotate that coordinate system
	// back to G4, not take vector in G4 and rotate to new position.  hence inverse
            pos = primaryProtonGunRotation.inverse()*pos + primaryProtonGunOrigin;
            momentum = primaryProtonGunRotation.inverse()*momentum;
      } else {
        edm::LogError("KINEMATICS")
          << "Do not know what to do with this generator id: " 
          << genpart._generatorId
          << "  Skipping this track.";
        continue;
      }

      // Shorthand for readability: 4-momentum.
	 HepLorentzVector p4(genpart._momentum);
	 p4.setV(momentum);
    
      // Create a new vertex 
      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos,genpart._time);
    
      // Create a particle.
      G4PrimaryParticle* particle = 
        new G4PrimaryParticle(genpart._pdgId, p4.x(), p4.y(), p4.z(), p4.e() );
    
      // Set the charge.  Do I really need to do this?
      G4ParticleDefinition const* g4id = particle->GetG4code();
      particle->SetCharge( g4id->GetPDGCharge()*eplus );
          
      // Add the particle to the event.
      vertex->SetPrimary( particle );
      
      // Add the vertex to the event.
      event->AddPrimaryVertex( vertex );

    }

  }


  // A very simple track for debugging G4 volumes and graphics
  void PrimaryGeneratorAction::testTrack(G4Event* event){ 

    // All tracks start from the same spot.
    G4ThreeVector const& position = _world->getMu2eDetectorOrigin();

    // Magnitude of the momentum.
    G4double p0  = 50. + 100.*G4UniformRand();

    // Generate the momentum 3-vector.
    G4ThreeVector momentum = p0 * _randomUnitSphere->shoot();

    /*
    // Status report.
    printf ( "Forward: %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f\n",
    1,
    position.x(),
    position.y(),
    position.z(),
    momentum.x(),
    momentum.y(),
    momentum.z(),
    p0
    );
    */
 
    // Create a new vertex 
    G4PrimaryVertex* vertex = new G4PrimaryVertex(position,0.);

    // Create a new particle.
    G4PrimaryParticle* particle = new 
      G4PrimaryParticle(_particleDefinition,
                        momentum.x(),
                        momentum.y(),
                        momentum.z()
                        );
    particle->SetMass( 0. );
    particle->SetCharge( eplus );
    particle->SetPolarization( 0., 0., 0.);

    // Add the particle to the event.
    vertex->SetPrimary( particle );
  
    // Add the vertex to the event.
    event->AddPrimaryVertex( vertex );

  }

}
