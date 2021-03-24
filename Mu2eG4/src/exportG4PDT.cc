//
// Export the G4 particle data table.
//
//
// Contact person Rob Kutschke
//

#include "Mu2eG4/inc/exportG4PDT.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "Geant4/G4ParticleTable.hh"

#include <iostream>

using namespace std;

namespace mu2e{

  void exportG4PDT( string const& tag ){

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();

    int nG4Extra(0);

    iter->reset();
    while( (*iter)() ){

      G4ParticleDefinition* particle = iter->value();
      int pdgId                      = particle->GetPDGEncoding ();
      G4String particleName          = particle->GetParticleName();
      cout << "exportG4PDT "
           << tag   << " "
           << pdgId << " "
           << particleName << " "
           << particle->GetPDGCharge() << "  "
           << particle->GetPDGMass() << "  "
           << particle->GetPDGWidth() << "  "
           << particle->GetParticleType() << "  "
           << particle->GetParticleSubType() << "  "
           << particle->GetPDGStable() << "  "
           << endl;
      if ( pdgId > PDGCode::G4Threshold ){
        ++nG4Extra;
      }

    }
    cout << "exportG4PDT "
         << tag << " "
         << " sizes: "
         << ptable->size() << " "
         << nG4Extra << endl;

  }

}  // end namespace mu2e
