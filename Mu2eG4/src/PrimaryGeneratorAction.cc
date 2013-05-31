//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.cc,v 1.46 2013/05/31 18:06:28 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:06:28 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <stdexcept>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// G4 Includes
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "globals.hh"

// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "MCDataProducts/inc/PDGCode.hh"

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

    // those should really be references; but because of the "if" they are not...

    G4ThreeVector mu2eOrigin;

    SimpleConfig const & _config = (*(art::ServiceHandle<GeometryService>())).config();

    // check if this is standard mu2e configuration; not all generators may work if it is not
    if (_config.getBool("mu2e.standardDetector",true)) {

      GeomHandle<WorldG4>  worldGeom;
      // Get the offsets to map from generator world to G4 world.
      mu2eOrigin      = worldGeom->mu2eOriginInWorld();
    }

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
      pos += mu2eOrigin;

      // Create a new vertex
      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos,genpart.time());

      G4int pPdgId = genpart.pdgId();

      static int const verbosityLevel = 0; // local debugging variable
      if ( verbosityLevel > 0) {
        cout << __func__ << " genpart.pdgId()   : " <<pPdgId << endl;
      }

//       static bool firstTime = true; // uncomment generate all nuclei ground states
//       if (firstTime) {
//         G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIon();
//         firstTime = false;
//       }

      if (pPdgId>PDGCode::G4Threshold) {

        G4int ZZ,AA,LL,JJ; 
        G4double EE;

        int exc = pPdgId%10;

        // subtract exc to get Z,A,...

        bool retCode = G4IonTable::GetNucleusByEncoding(pPdgId-exc,ZZ,AA,LL,EE,JJ);

        // if g4 complains about no GenericIon we need to abort and add it in the physics list...

        if ( verbosityLevel > 0) {
          cout << __func__ << " will set up nucleus pPdgId,Z,A,L,E,J, exc: " 
               << pPdgId 
               << ", " << ZZ << ", " << AA << ", " << LL << ", " << EE << ", " << JJ
               << ", " << exc << ", " << retCode
               << endl;
        } 

        G4ParticleDefinition const * pDef;

        if (exc==0) {

          // a ground state

          pDef = G4ParticleTable::GetParticleTable()->GetIon(ZZ,AA,LL,0.0);
          // looks like spin is ignored and should never be non zero...

        } else {

          // an excited state; we will fudge it by adding 1keV
          pDef = G4ParticleTable::GetParticleTable()->GetIon(ZZ,AA,LL,0.001);

        }

        if ( verbosityLevel > 0) {
          cout << __func__ << " will set up : " 
               << pDef->GetParticleName() 
               << " with id: " << pDef->GetPDGEncoding() 
               << endl;
        }
        // the ids should be the same
        if ( pPdgId != pDef->GetPDGEncoding() ) {

          throw cet::exception("GENE")
            << "Problem creating "<<pPdgId << "\n";
        }

      }

      G4PrimaryParticle* particle =
        new G4PrimaryParticle(pPdgId,
                              momentum.x(),
                              momentum.y(),
                              momentum.z(),
                              genpart.momentum().e() );

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
    GeomHandle<WorldG4>  world;
    GeomHandle<Mu2eBuilding>  building;
    G4ThreeVector const& position = building->relicMECOOriginInMu2e() + world->mu2eOriginInWorld();

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
