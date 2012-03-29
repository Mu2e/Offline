//
// Export the G4 particle data table.
//
// $Id: exportG4PDT.cc,v 1.1 2012/03/29 17:03:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/03/29 17:03:14 $
//
// Contact person Rob Kutschke
//

#include "Mu2eG4/inc/exportG4PDT.hh"
#include "MCDataProducts/inc/PDGCode.hh"

#include "G4ParticleTable.hh"

#include <iostream>

using namespace std;

namespace mu2e{

  void exportG4PDT( ){

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();

    int nG4Extra(0);

    cout << "exportG4PDT: " << endl;
    iter->reset();
    while( (*iter)() ){

      G4ParticleDefinition* particle = iter->value();
      int pdgId                      = particle->GetPDGEncoding ();
      G4String particleName          = particle->GetParticleName();
      cout << "exportG4PDT: "
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
    cout << "exportG4PDT: sizes: " << ptable->size() << " " << nG4Extra << endl;

  }

}  // end namespace mu2e
