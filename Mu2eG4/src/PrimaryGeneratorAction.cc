//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.cc,v 1.51 2013/10/05 05:04:48 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/05 05:04:48 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <stdexcept>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Utilities/Globals.h"

// G4 Includes
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#if G4VERSION<4099
#include "G4ParticleTable.hh"
#endif
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4Threading.hh"

// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"


using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

namespace mu2e {

  PrimaryGeneratorAction::PrimaryGeneratorAction(const Mu2eG4Config::Debug& debug,
                                                 Mu2eG4PerThreadStorage* tls)
    : verbosityLevel_(debug.diagLevel())
    , perThreadObjects_(tls)
  {
    if ( verbosityLevel_ > 0 ) {
      cout << __func__ << " verbosityLevel_  : " <<  verbosityLevel_ << endl;
    }
  }

  //load in per-art-event data from GenEventBroker and per-G4-event data from EventObjectManager
  void PrimaryGeneratorAction::setEventData()
  {
    genParticles_ = perThreadObjects_->gensHandle.isValid() ?
      perThreadObjects_->gensHandle.product() :
      nullptr;

    hitInputs_ = perThreadObjects_->genInputHits;

    parentMapping_ = perThreadObjects_->simParticlePrimaryHelper;

  }


  void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
  {
    //load the GenParticleCollection etc for the event
    setEventData();
    // For debugging.
    //testTrack(anEvent);

    fromEvent(anEvent);
  }

  // Copy generated particles from the event into G4.
  // At the moment none of the supported generators make multi-particle vertices.
  // That may change in the future.
  void PrimaryGeneratorAction::fromEvent(G4Event* event)
  {

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const & _config = geom->config();

    // Get the offsets to map from generator world to G4 world.
    // check if this is standard mu2e configuration or not

    G4ThreeVector const mu2eOrigin =
      (!_config.getBool("mu2e.standardDetector",true) || !(geom->isStandardMu2eDetector()))
      ?  G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();

    // For each generated particle, add it to the event.
    if(genParticles_) {
      for (unsigned i=0; i < genParticles_->size(); ++i) {
        const GenParticle& genpart = (*genParticles_)[i];
        addG4Particle(event,
                      genpart.pdgId(),
                      // Transform into G4 world coordinate system
                      genpart.position() + mu2eOrigin,
                      genpart.time(),
                      genpart.properTime(),
                      genpart.momentum());

        parentMapping_->addEntryFromGenParticle(i);
      }
    }


    // a test
    // int ovl = verbosityLevel_;
    // verbosityLevel_ = 2;
    // addG4Particle(event,
    //               static_cast<PDGCode::type>(1000010048),
    //               // Transform into G4 world coordinate system
    //               G4ThreeVector()+mu2eOrigin,
    //               0.0,
    //               0.0,
    //               G4ThreeVector());
    // verbosityLevel_ = ovl;


    // Also create particles from the input hits
    for(const auto& hitcoll : *hitInputs_) {

      for(const auto& hit : *hitcoll) {
        addG4Particle(event,
                      hit.simParticle()->pdgId(),
                      // Transform into G4 world coordinate system
                      hit.position() + mu2eOrigin,
                      hit.time(),
                      hit.properTime(),
                      hit.momentum());

        parentMapping_->addEntryFromStepPointMC(hit.simParticle()->id());

      }
    }

  }


  void PrimaryGeneratorAction::addG4Particle(G4Event *event,
                                             PDGCode::type pdgId,
                                             const G4ThreeVector& pos,
                                             double time,
                                             double properTime,
                                             const G4ThreeVector& mom)
  {
    // Create a new vertex
    G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

    if ( verbosityLevel_ > 1) {
      cout << __func__ << " pdgId   : " <<pdgId << endl;
    }

    //       static bool firstTime = true; // uncomment generate all nuclei ground states
    //       if (firstTime) {
    //#if G4VERSION<4099
    //         G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIon();
    //#else
    //         G4IonTable::GetIonTable()->GetIonTable()->CreateAllIon();
    //#endif
    //         firstTime = false;
    //       }

    if (pdgId>PDGCode::G4Threshold) {

      G4int ZZ,AA,LL,JJ;
      G4double EE;

      int exc = pdgId%10;

      // subtract exc to get Z,A,...

      bool retCode = G4IonTable::GetNucleusByEncoding(pdgId-exc,ZZ,AA,LL,EE,JJ);

      // if g4 complains about no GenericIon we need to abort and add it in the physics list...

      if ( verbosityLevel_ > 1) {
        cout << __func__ << " will set up nucleus pdgId,Z,A,L,E,J, exc: "
             << pdgId
             << ", " << ZZ << ", " << AA << ", " << LL << ", " << EE << ", " << JJ
             << ", " << exc << ", " << retCode
             << endl;
      }

      G4ParticleDefinition const * pDef;

      G4double excEnergy = 0.001;
      // for an excited state; we will fudge it by adding 1keV regardles of the isomer level

      if (exc==0) {

        // a ground state

#if G4VERSION<4099
        pDef = G4ParticleTable::GetParticleTable()->GetIon(ZZ,AA,LL,0.0);
#else
        pDef = G4IonTable::GetIonTable()->GetIon(ZZ,AA,LL,0.0);
#endif
        // looks like spin is ignored and should never be non zero...

      } else {

#if G4VERSION<4099
        pDef = G4ParticleTable::GetParticleTable()->GetIon(ZZ,AA,LL,excEnergy);
#else
        pDef = G4IonTable::GetIonTable()->GetIon(ZZ,AA,LL,excEnergy);
#endif
      }

      if ( verbosityLevel_ > 1) {
        cout << __func__ << " will set up : "
             << pDef->GetParticleName()
             << " with id: " << pDef->GetPDGEncoding()
             << endl;
      }
      // the ids should be the same
      if ( pdgId != pDef->GetPDGEncoding() ) {

        throw cet::exception("GENE")
          << "Problem creating " << pdgId << "\n";
      }

    }

    G4PrimaryParticle* particle =
      new G4PrimaryParticle(pdgId,
                            mom.x(),
                            mom.y(),
                            mom.z());

    // Add the particle to the event.
    vertex->SetPrimary( particle );

    // Add the vertex to the event.
    event->AddPrimaryVertex( vertex );
  }


} // end namespace mu2e
