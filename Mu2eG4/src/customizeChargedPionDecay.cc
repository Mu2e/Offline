//
// Add the decay pi+ -> e+ nu_e to the G4 decay table for pi+ and similarly
// for pi-.
//
// $Id: customizeChargedPionDecay.cc,v 1.1 2012/07/10 21:16:53 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/10 21:16:53 $
//
// In default configured G4, the only decy mode for pi+ is pi+ -> mu+ nu_mu.
// Similarly for pi-.  This function looks at the config file parameter
// g4.PiENuPolicy, and takes one of the following actions, depending on
// the value of that parameter:
//
//  Value    Action
//   PDG   - Set the pi -> e nu and pi -> mu nu branching fractions to
//           their PDG values.
//   None  - Leave the decay table unchanged ( 100% pi-> mu nu ).
//   All   - Set the pi-> e nu branching fraction to 100%
//   nnnnn - where nnn is a numerical value in the range [0., 1.]
//           Set the e nu branching fraction to the specified number.
//           and the mu nu branching fraction to (1-nnnnn).
//
// If the parameter g4.PiENuPolicy is absent, the default is PDG.
//
// In all cases set B( pi+ -> mu+ nu_mu ) is 1. - B( pi+ -> e+ nu_e)
//
// This routine always modifies the decay tables for both pi+ and pi-.
//
// There is no way to remove the pre-existing pi+ -> mu+ nu_mu channel but
// we can set its BR to 0.
//
// The change to the mu nu BR must be made before the e nu BR is added.
// Parts of G4 require that the decay table be in order of increasing BR
// but the list is not resorted if a BR is changed.
//

#include "Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "cetlib/exception.h"

#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include <sstream>
#include <algorithm>

using namespace std;

namespace mu2e{

  namespace {

    // A helper class to describe one decay channel.
    struct Channel{

      Channel ( G4ParticleTable *pdt,
                PDGCode::type aparent,
                PDGCode::type achild0,
                PDGCode::type achild1
                ):
        parent(aparent),
        child0(achild0),
        child1(achild1),
        parentName(pdt->FindParticle(parent)->GetParticleName()),
        child0Name(pdt->FindParticle(child0)->GetParticleName()),
        child1Name(pdt->FindParticle(child1)->GetParticleName())
      {}

      PDGCode::type parent;
      PDGCode::type child0;
      PDGCode::type child1;
      G4String      parentName;
      G4String      child0Name;
      G4String      child1Name;

    }; // end class Channel

    // A helper class to parse the configuration parameter.
    struct Control {

      Control ( const SimpleConfig& config):
        brENu(0.),
        brMuNu(1.),
        addENu(true){

        // FIXME:  Get this from the global conditions service.
        double pdgBR = 1.23e-4;

        string name("g4.PiENuPolicy");

        string option = config.getString(name,"PDG");
        if ( option == "PDG" ){
          brENu = pdgBR;
          brMuNu = 1.-brENu;
          return;
        } else if ( option == "All" ){
          brENu  = 1.0;
          brMuNu = 0.0;
        } else if ( option == "None" ){
          addENu = false;
        } else {

          // The option must be interpretable as a number on the range [0,1] .
          istringstream is(option);
          double b;
          is >> b;

          // Not interpretable as a number.
          if ( !is ){
            throw cet::exception("G4")
              << "customizeChargedPionDecay: could not parse the option value for " << name << ": \n"
              << "    the option value was: " << option << "\n";
          }

          // Out of range.
          if ( b < 0. || b > 1. ){
            throw cet::exception("G4")
              << "customizeChargedPionDecay: branching ratio specified by " << name << " is out of range: \n"
              << "    branching ratio was: " << b << " (" << option << " )\n";
          }

          // Treat this the same as "None"
          if ( b == 0. ){
            addENu = false;
            return;
          }

          // Passed all of the checks, so accept the br from the config file.
          // Paranoia about round off errors: hence std::max
          brENu  = b;
          brMuNu = std::max(1.-brENu, 0.);

        }
      }

      // Branching fraction for pi -> e nu_e
      double brENu;

      // Branching fraction for pi -> mu nu_mu
      double brMuNu;

      // Do we need to add the pi -> e nu_e channel to the list of decay channels?
      bool addENu;

    }; // end class Control

  } // end anonymous namespace

  void customizeChargedPionDecay(const SimpleConfig& config) {

    // Figure out what we need to do.
    Control ctrl(config);
    int verbosity = config.getInt("g4.PiENuPolicyVerbosity",0);

    // Nothing to do or to print out.
    if ( !ctrl.addENu && verbosity == 0 ) return;

    G4ParticleTable *pdt = G4ParticleTable::GetParticleTable();

    // Define the decay channels we will create.
    std::vector<Channel> pions;
    pions.push_back( Channel( pdt, PDGCode::pi_plus,  PDGCode::e_plus,  PDGCode::nu_e      ) );
    pions.push_back( Channel( pdt, PDGCode::pi_minus, PDGCode::e_minus, PDGCode::anti_nu_e ) );

    for ( std::vector<Channel>::const_iterator i=pions.begin();
          i != pions.end(); ++i ) {

      PDGCode::type parent            = i->parent;
      G4ParticleDefinition * particle = pdt->FindParticle( parent );
      G4DecayTable * decayTable       = particle->GetDecayTable();

      if ( decayTable->entries() != 1 ){
        throw cet::exception("G4")
          << "customizeChargedPionDecay: expected to find only one entry in the decay table for " << parent
          << " (" << particle->GetParticleName() << ")\n"
          << "Number found: " << decayTable->entries() << "\n";
      }

      if ( ctrl.addENu ) {

        // Modify BR of the mu nu_mu mode.
        G4VDecayChannel* muNuChannel = (*decayTable)[0];
        muNuChannel->SetBR(ctrl.brMuNu);

        // Add the e nu_e mode
        G4VDecayChannel* mode =
          new G4PhaseSpaceDecayChannel( i->parentName, ctrl.brENu, 2, i->child0Name, i->child1Name );
        decayTable->Insert(mode);
      }

      if ( verbosity > 0 ) {
        decayTable->DumpInfo();
      }

    }

  }

} // end namespace mu2e
