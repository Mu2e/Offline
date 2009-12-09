//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.cc,v 1.3 2009/12/09 18:56:10 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2009/12/09 18:56:10 $
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
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#include <boost/shared_ptr.hpp>

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

PrimaryGeneratorAction::PrimaryGeneratorAction()
  :_randomUnitSphere(0)
{
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
// At the momment none of the supported generators make multi-particle vertices.
// That may change in the future.
// 
void PrimaryGeneratorAction::fromEvent(G4Event* event){ 

  // Get the offsets to map from generator world to G4 world.
  G4ThreeVector const& cosmicReferencePlane = _world->getCosmicReferencePoint();
  G4ThreeVector const& detectorOrigin       = _world->getMu2eDetectorOrigin();

  // Get generated particles from the event.
  edm::Handle<ToyGenParticleCollection> handle;
  _event->getByLabel("generate",handle);

  // Fill multiplicity histogram.
  _totalMultiplicity->Fill( handle->size() );
  
  // For each generated particle, add it to the event.
  for ( unsigned int i=0; i<handle->size(); ++i){

    ToyGenParticle const& genpart = (*handle)[i];

    // Adjust the origin into the G4 world.
    G4ThreeVector     pos(genpart._position);
    if( genpart._generatorId == GenId::conversionGun ||
	genpart._generatorId == GenId::dio1 ||
	genpart._generatorId == GenId::pionCapture  ){
      	pos += detectorOrigin;
    } else if ( genpart._generatorId == GenId::cosmicToy ){
      pos += cosmicReferencePlane;
    } else {
      edm::LogError("KINEMATICS")
	<< "Do not know what to do with this generator id: " 
	<< genpart._generatorId
	<< "  Skipping this track.";
      continue;
    }

    // Shorthand for readability: 4-momentum.
    HepLorentzVector const& p4(genpart._momentum);
    
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


// A very simple track for debugging G4 volumes and graphics.
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
