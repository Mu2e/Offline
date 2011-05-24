//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.cc,v 1.23 2011/05/24 17:19:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:19:03 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
//#include <cassert>
#include <stdexcept>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// G4 Includes
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "globals.hh"

// Mu2e includes
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// ROOT includes
#include "TH1D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

  PrimaryGeneratorAction::PrimaryGeneratorAction( const string& generatorModuleLabel ):
    _generatorModuleLabel(generatorModuleLabel){

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    _particleDefinition = particleTable->FindParticle("chargedgeantino");

    // Book histograms.
    art::ServiceHandle<art::TFileService> tfs;
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
    G4ThreeVector const& mu2eOrigin                  = _world->getMu2eOrigin();
    G4ThreeVector const& cosmicReferencePlane        = _world->getCosmicReferencePoint();
    G4ThreeVector const& detectorOrigin              = _world->getMu2eDetectorOrigin();
    G4ThreeVector const& primaryProtonGunOrigin      = _world->getPrimaryProtonGunOrigin();
    G4RotationMatrix const& primaryProtonGunRotation = _world->getPrimaryProtonGunRotation();

    // Get generated particles from the event.
    art::Handle<GenParticleCollection> handle;
    _event->getByLabel(_generatorModuleLabel,handle);

    // Fill multiplicity histogram.
    _totalMultiplicity->Fill( handle->size() );

    // For each generated particle, add it to the event.
    for ( unsigned int i=0; i<handle->size(); ++i){

      GenParticle const& genpart = (*handle)[i];

      // Transform from generator coordinate system G4 world coordinate system.
      G4ThreeVector      pos(genpart.position());
      G4ThreeVector momentum(genpart.momentum().v());

      if( genpart.generatorId() == GenId::conversionGun     ||
          genpart.generatorId() == GenId::dio1              ||
          genpart.generatorId() == GenId::ejectedProtonGun  ||
          genpart.generatorId() == GenId::ejectedNeutronGun ||
          genpart.generatorId() == GenId::pionCapture       ||
          genpart.generatorId() == GenId::piEplusNuGun      ||
          genpart.generatorId() == GenId::nuclearCaptureGun ){
        pos += detectorOrigin;
      } else if ( genpart.generatorId() == GenId::cosmicToy ||
                  genpart.generatorId() == GenId::cosmicDYB ||
                  genpart.generatorId() == GenId::cosmic ){
        pos += cosmicReferencePlane;
      } else if ( genpart.generatorId() == GenId::primaryProtonGun ){
        pos = primaryProtonGunRotation*pos + primaryProtonGunOrigin;
        momentum = primaryProtonGunRotation*momentum;
      } else if ( genpart.generatorId() == GenId::particleGun ||
                  genpart.generatorId() == GenId::fromG4BLFile ){
        pos += mu2eOrigin;
      } else {
        mf::LogError("KINEMATICS")
          << "Do not know what to do with this generator id: "
          << genpart.generatorId()
          << "  Skipping this track.";
        continue;
      }

      // Create a new vertex
      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos,genpart.time());

      G4PrimaryParticle* particle =
        new G4PrimaryParticle(genpart.pdgId(),
                              momentum.x(),
                              momentum.y(),
                              momentum.z(),
                              genpart.momentum().e() );

      // Set the charge.  Do I really need to do this?
      G4ParticleDefinition const* g4id = particle->GetG4code();
      particle->SetCharge( g4id->GetPDGCharge()*eplus );

      // Add the particle to the event.
      vertex->SetPrimary( particle );

      // Add the vertex to the event.
      event->AddPrimaryVertex( vertex );

    }

  }

  // A very simple generator for debugging G4 volumes and graphics.
  void PrimaryGeneratorAction::testTrack(G4Event* event){

    // Generator for uniform coverage of a restricted region on a unit sphere.
    static RandomUnitSphere randomUnitSphere( *CLHEP::HepRandom::getTheEngine(), -0.7, 0.7 );

    // All tracks start from the same spot.
    G4ThreeVector const& position = _world->getMu2eDetectorOrigin();

    // Magnitude of the momentum.
    G4double p0  = 50. + 100.*G4UniformRand();

    // Generate the momentum 3-vector.
    G4ThreeVector momentum = randomUnitSphere.fire(p0);

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

} // end namespace mu2e
