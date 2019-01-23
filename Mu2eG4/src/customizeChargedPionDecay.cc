//
// Add the decay pi+ -> e+ nu_e to the G4 decay table for pi+ and similarly
// for pi-.
//
// $Id: customizeChargedPionDecay.cc,v 1.3 2012/07/20 00:26:53 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/20 00:26:53 $
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
// Notes:
// 1) If the parameter g4.PiENuPolicy is absent, the default is PDG.
//
// 2) In all cases set B( pi+ -> mu+ nu_mu ) is 1. - B( pi+ -> e+ nu_e)
//
// 3) This routine always modifies the decay tables for both pi+ and pi-.
//
// 4) There is no way to remove the pre-existing pi+ -> mu+ nu_mu channel but
//    we can set its BR to 0.
//
// 5) The change to the mu nu BR must be made before the e nu BR is added.
//    Parts of G4 require that the decay table be in order of increasing BR
//    and an addition causes a resort. But changing a BR does not, itself,
//    trigger a resort.
//
// 6) If we are running with a less than complete physics list, for example
//    transportOnly then some particle or particle table information will
//    be incomplete.  Test for those cases and skip them.
//

#include "Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include <sstream>
#include <algorithm>

#include <iostream>

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
        ok(true),
        parent(aparent),
        child0(achild0),
        child1(achild1),
        parentName(name(pdt,parent)),
        child0Name(name(pdt,child0)),
        child1Name(name(pdt,child1)){
      }

      // Is all particle information present; see note 6.
      bool ok;

      // Decay channel information
      PDGCode::type parent;
      PDGCode::type child0;
      PDGCode::type child1;
      G4String      parentName;
      G4String      child0Name;
      G4String      child1Name;

      // Get name string for the particle. See note 6.
      G4String name( G4ParticleTable * pdt, PDGCode::type pid ){
        G4ParticleDefinition * particle  = pdt->FindParticle(pid);
        if ( particle == 0 ){
          ok = false;
          mf::LogWarning("G4") << "There is no G4 particle information for PDGcode: " << pid
                               << "\nHope that's OK. Skipping and continuing ...\n";
          return G4String();
        }
        return particle->GetParticleName();
      }

    }; // end class Channel

    // A helper class to parse the configuration parameter.
    struct Control {

      Control (const std::string& option):
        brENu(0.),
        brMuNu(1.),
        addENu(true){

        // FIXME:  Get this from the global conditions service.
        double pdgBR = 1.23e-4;

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
              << "customizeChargedPionDecay: could not parse the option value = "
              << option << "\n";
          }

          // Out of range.
          if ( b < 0. || b > 1. ){
            throw cet::exception("G4")
              << "customizeChargedPionDecay: branching ratio is out of range: \n"
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

      Control ( const SimpleConfig& config)
        : Control(config.getString("g4.PiENuPolicy","PDG"))
      {}

      Control ( const fhicl::ParameterSet& pset)
        : Control(pset.get<std::string>("physics.PiENuPolicy"))
      {}

      // Branching fraction for pi -> e nu_e
      double brENu;

      // Branching fraction for pi -> mu nu_mu
      double brMuNu;

      // Do we need to add the pi -> e nu_e channel to the list of decay channels?
      bool addENu;

    }; // end class Control

    int getPiENuPolicyVerbosity(const SimpleConfig& config) {
      return config.getInt("g4.PiENuPolicyVerbosity",0);
    }
    int getPiENuPolicyVerbosity(const fhicl::ParameterSet& pset) {
      return pset.get<int>("debug.PiENuPolicyVerbosity");
    }

  } // end anonymous namespace

  template<class Config>
  void customizeChargedPionDecay(const Config& config) {

    // Figure out what we need to do.
    Control ctrl(config);
    int verbosity = getPiENuPolicyVerbosity(config);

    // Nothing to do or to print out.
    if ( !ctrl.addENu && verbosity == 0 ) return;

    G4ParticleTable *pdt = G4ParticleTable::GetParticleTable();

    // Define the decay channels we will create.
    std::vector<Channel> pions;
    pions.push_back( Channel( pdt, PDGCode::pi_plus,  PDGCode::e_plus,  PDGCode::nu_e      ) );
    pions.push_back( Channel( pdt, PDGCode::pi_minus, PDGCode::e_minus, PDGCode::anti_nu_e ) );

    for ( std::vector<Channel>::const_iterator i=pions.begin();
          i != pions.end(); ++i ) {

      // Incomplete information - skip this channel. See note 6.
      if ( ! i-> ok ) return;

      PDGCode::type parent            = i->parent;
      G4ParticleDefinition * particle = pdt->FindParticle( parent );

      // See note 6.
      if ( particle  == 0 ){
        mf::LogWarning("G4") << "There is no G4 particle information for PDGcode: " << parent
                             << "\nHope that's OK. Skipping and continuing ...\n";
        continue;
      }

      G4DecayTable * decayTable = particle->GetDecayTable();

      // See note 6.
      if ( decayTable == 0 ){
        mf::LogWarning("G4") << "There is no decay table for PDGcode: " << parent
                             << "\nHope that's OK. Skipping and continuing ...\n";
        continue;
      }


      if ( decayTable->entries() != 1 ){

        G4cout
          << __func__ << " : expected to find only one entry in the decay table for " << parent
          << " (" << particle->GetParticleName() << ")\n"
          << "Number found: " << decayTable->entries() << "\n";
        decayTable->DumpInfo();

        throw cet::exception("G4")
          << __func__ << " : expected to find only one entry in the decay table for " << parent
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

  template void customizeChargedPionDecay(const SimpleConfig&);
  template void customizeChargedPionDecay(const fhicl::ParameterSet&);

} // end namespace mu2e
