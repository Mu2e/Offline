//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
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
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4Threading.hh"
#include "G4IsotopeProperty.hh"
#include "G4NuclideTable.hh"

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

namespace mu2e {

  PrimaryGeneratorAction::PrimaryGeneratorAction(const Mu2eG4Config::Debug& debug,
                                                 Mu2eG4PerThreadStorage* tls)
    : verbosityLevel_(debug.diagLevel())
    , testPDGIdToGenerate_(debug.ionToGenerate(testIonToGenerate_))
    , perThreadObjects_(tls)
  {
    if ( verbosityLevel_ > 0 ) {
      G4cout << __func__ << " verbosityLevel_  : " <<  verbosityLevel_ << G4endl;
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
                      0.0, // Mu2e GenParticles are not excited ions
                      0,
                      // Transform into G4 world coordinate system
                      genpart.position() + mu2eOrigin,
                      genpart.time(),
                      genpart.properTime(),
                      genpart.momentum());

        parentMapping_->addEntryFromGenParticle(i);
      }
    }

    if ( !testPDGIdToGenerate_ ) {

      // standard flow

      // Also create particles from the input hits
      for(const auto& hitcoll : *hitInputs_) {

        for(const auto& hit : *hitcoll) {
          addG4Particle(event,
                        hit.simParticle()->pdgId(),
                        hit.simParticle()->startExcitationEnergy(),
                        hit.simParticle()->startFloatLevelBaseIndex(),
                        // Transform into G4 world coordinate system
                        hit.position() + mu2eOrigin,
                        hit.time(),
                        hit.properTime(),
                        hit.momentum());

          parentMapping_->addEntryFromStepPointMC(hit.simParticle()->id());

        }
      }

    } else {

      // testing isomer creation; will fail shortly after the
      // creation, but will print a lot of diagnostic info
      PDGCode::type iid = static_cast<PDGCode::type>(std::get<0>(testIonToGenerate_));
      G4cout << __func__ << " *TESTING*: Adding *ONLY*  : " << iid << G4endl;

      // G4IonTable::GetIonTable()->CreateAllIsomer();
      // the above does not create them

      // G4NuclideTable::GetNuclideTable()->GenerateNuclide();

      // G4cout << __func__ << " Nuclide table " << G4endl;
      // G4NuclideTable::GetNuclideTable()->DumpTable();

      // G4cout << __func__ << " Ion table " << G4endl;
      // G4IonTable::GetIonTable()->DumpTable("ALL");

      int ovl = verbosityLevel_;
      int govl = G4ParticleTable::GetParticleTable()->GetVerboseLevel();
      verbosityLevel_ = 2;
      G4ParticleTable::GetParticleTable()->SetVerboseLevel(verbosityLevel_);
      addG4Particle(event,
                    iid,
                    (iid%10>0) ? std::get<1>(testIonToGenerate_) : 0.0,
                    (iid%10>0) ? std::get<2>(testIonToGenerate_) : 0,
                    // Transform into G4 world coordinate system
                    G4ThreeVector()+mu2eOrigin,
                    0.0,
                    0.0,
                    G4ThreeVector());
      verbosityLevel_ = ovl;
      G4ParticleTable::GetParticleTable()->SetVerboseLevel(govl);

    }
  }

  void PrimaryGeneratorAction::addG4Particle(G4Event *event,
                                             PDGCode::type pdgId,
                                             double excitationEnergy,
                                             int floatLevelBaseIndex,
                                             const G4ThreeVector& pos,
                                             double time,
                                             double properTime,
                                             const G4ThreeVector& mom)
  {
    // Create a new vertex
    G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

    if ( verbosityLevel_ > 1) {
      G4cout << __func__ << " Setting up pdgId   : " <<pdgId << G4endl;
    }

    G4ParticleDefinition const * pDef(nullptr);

    if (pdgId>PDGCode::G4Threshold) {

      // an ion

      G4int ZZ,AA, lvl;
      G4double EE(0.0);
      // get A & A based on pdgId
      bool retCode = G4IonTable::GetNucleusByEncoding(pdgId,ZZ,AA,EE,lvl);

      // if g4 complains about no GenericIon we need to abort and add it in the physics list...

      if ( verbosityLevel_ > 1) {
        G4cout << __func__ << " preparing nucleus pdgId,Z,A,lvl,ret: "
               << pdgId
               << ", " << ZZ << ", " << AA
               << ", " << lvl << ", " << retCode
               << G4endl;
      }

      pDef = G4IonTable::GetIonTable()->GetIon(ZZ, AA, excitationEnergy,
                                               G4Ions::FloatLevelBase(floatLevelBaseIndex));

      if (verbosityLevel_ > 1) {
        pDef->DumpTable();
      }

      // the ids should be the same
      if ( !pDef || pdgId != pDef->GetPDGEncoding() ) {
        if (pDef) {G4cout << __func__ << " got " << pDef->GetPDGEncoding() << " expected " << pdgId << G4endl;}
        throw cet::exception("GENE")
          << "Problem creating " << pdgId << "\n";

      }

      if ( verbosityLevel_ > 1) {
        G4cout << __func__ << " will set up : "
               << pDef->GetParticleName()
               << " with id: " << pDef->GetPDGEncoding()
               << G4endl;
      }

    }

    // note that the particle definition not the pdg code is used for ions
    // Add the particle to the event.
    vertex->SetPrimary( pDef ?
                        new G4PrimaryParticle(pDef, mom.x(),mom.y(),mom.z()) :
                        new G4PrimaryParticle(pdgId,mom.x(),mom.y(),mom.z()) );

    // Add the vertex to the event.
    event->AddPrimaryVertex( vertex );
  }


} // end namespace mu2e
