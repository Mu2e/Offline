//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Utilities/Globals.h"

// G4 Includes
#include "Geant4/G4Event.hh"
#include "Geant4/G4IonTable.hh"
#include "Geant4/globals.hh"
#include "Geant4/G4Threading.hh"
#include "Geant4/G4IsotopeProperty.hh"
#include "Geant4/G4NuclideTable.hh"

// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PrimaryGeneratorAction.hh"
#include "Offline/Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "Offline/Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Offline/Mu2eG4/inc/SimParticleHelper.hh"


using namespace std;

namespace mu2e {

  Mu2eG4PrimaryGeneratorAction::Mu2eG4PrimaryGeneratorAction(const Mu2eG4Config::Debug& debug,
                                                 Mu2eG4PerThreadStorage* tls)
    : verbosityLevel_(debug.diagLevel())
    , testPDGIdToGenerate_(debug.ionToGenerate(testIonToGenerate_))
    , perThreadObjects_(tls)
  {
    if ( verbosityLevel_ > 0 ) {
      G4cout << __func__ << " verbosityLevel_  : " <<  verbosityLevel_ << G4endl;
    }
  }


  void Mu2eG4PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const & _config = geom->config();

    // Get the offsets to map from generator world to G4 world.
    // check if this is standard mu2e configuration or not

    // FIXME: is the SimpleConfig access here thread safe?
    G4ThreeVector const mu2eOrigin =
      (!_config.getBool("mu2e.standardDetector",true) || !(geom->isStandardMu2eDetector()))
      ?  G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();


    auto artEvent = perThreadObjects_->artEvent;
    auto& inputs = perThreadObjects_->ioconf.inputs();

    switch(inputs.primaryType().id()) {
    default: throw cet::exception("CONFIG")
        << "Error: Mu2eG4PrimaryGeneratorAction: unknown Mu2eG4 primaryType id = "
        <<inputs.primaryType().id()
        <<std::endl;

    case Mu2eG4PrimaryType::GenParticles: {
      auto const h = artEvent->getValidHandle<GenParticleCollection>(inputs.primaryTag());
      for (unsigned i=0; i < h->size(); ++i) {
        const GenParticle& genpart = (*h)[i];
        addG4Particle(event,
                      genpart.pdgId(),
                      0.0, // Mu2e GenParticles are not excited ions
                      0,
                      // Transform into G4 world coordinate system
                      genpart.position() + mu2eOrigin,
                      genpart.time(),
                      genpart.properTime(),
                      genpart.momentum());

        perThreadObjects_->simParticlePrimaryHelper->addEntry(art::Ptr<GenParticle>(h, i));
      }
    }
      break; // GenParticles

    case Mu2eG4PrimaryType::StepPoints: {
      auto const h = artEvent->getValidHandle<StepPointMCCollection>(inputs.primaryTag());
      for(const auto& hit : *h) {
        addG4Particle(event,
                      hit.simParticle()->pdgId(),
                      hit.simParticle()->startExcitationEnergy(),
                      hit.simParticle()->startFloatLevelBaseIndex(),
                      // Transform into G4 world coordinate system
                      hit.position() + mu2eOrigin,
                      hit.time(),
                      hit.properTime(),
                      hit.momentum());

        perThreadObjects_->simParticlePrimaryHelper->addEntry(&hit);
      }
    }
      break; // StepPoints

    case Mu2eG4PrimaryType::SimParticleLeaves: {
      auto const h = artEvent->getValidHandle<SimParticleCollection>(inputs.primaryTag());
      for(const auto& i : *h) {
        const SimParticle& particle = i.second;
        // Take the leaves only, we do not want to re-simulate
        // everything starting with the primary proton.
        if(particle.daughters().empty()) {
          addG4Particle(event,
                        particle.pdgId(),
                        particle.startExcitationEnergy(),
                        particle.startFloatLevelBaseIndex(),
                        // Transform into G4 world coordinate system
                        particle.endPosition() + mu2eOrigin,
                        particle.endGlobalTime(),
                        particle.endProperTime(),
                        particle.endMomentum());

          perThreadObjects_->simParticlePrimaryHelper->addEntry(&particle);
        }
      }
    }
      break; // SimParticles

    case Mu2eG4PrimaryType::StageParticles: {
      auto const h = artEvent->getValidHandle<StageParticleCollection>(inputs.primaryTag());
      for(const auto& s : *h) {
        addG4Particle(event,
                      s.pdgId(),

                      0.0, // no excited ions here
                      0,

                      // Transform into G4 world coordinate system
                      s.position() + mu2eOrigin,
                      s.time(),
                      0, //proper
                      s.momentum());

        perThreadObjects_->simParticlePrimaryHelper->addEntry(&s);
      }
    }
      break; // SimParticles
    }

    //----------------------------------------------------------------
    if ( testPDGIdToGenerate_ ) {

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
    //----------------------------------------------------------------

  }

  void Mu2eG4PrimaryGeneratorAction::addG4Particle(G4Event *event,
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
    else if (pdgId == PDGCode::geantino) {
      // AE: for some reason GetParticleTable()->FindParticle(0) does not find geantino
      // but GetParticleTable()->FindParticle("geantino") does
      pDef = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
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
