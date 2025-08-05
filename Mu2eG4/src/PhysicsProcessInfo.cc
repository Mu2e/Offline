//
// Information about physics processes.
//
//
// Original author Rob Kutschke
//
// Notes
// 1) For notes about the iterator class:
//      G4ParticleTable::G4PTblDicIterator
//    See http://mu2e.fnal.gov/atwork/computing/G4Notes.shtml .
//

// C++ includes
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <utility>

// Framework
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/PhysicsProcessInfo.hh"

// G4 includes
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4ProcessType.hh"

using namespace std;

namespace mu2e{


  PhysicsProcessInfo::PhysicsProcessInfo():
    _allProcesses(),
    _longestName(0){
  }

  void PhysicsProcessInfo::beginRun(){

    _allProcesses.clear();

    // Number of processes that are not known to the ProcessCode enum.
    int nUnknownProcesses(0);

    // a special process
    G4String gammaGP = G4String("GammaGeneralProc");

    // Get an iterator over existing particles. See note 1.
    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();
    iter->reset();

    // Loop over the known particles.
    while( (*iter)() ){

      // Find the list of processes attached to this particle.
      G4ParticleDefinition*  particle     = iter->value();
      G4ProcessManager*      pmanager     = particle->GetProcessManager();

      if ( !pmanager ){
        cout << __func__
             << " No process manager found for : " << particle->GetParticleName()
             << " skipping it "
             << endl;
        continue;
      }

      G4ProcessVector const* pVector      = pmanager->GetProcessList();
      G4String               particleName = particle->GetParticleName();

      // Loop over all processes.
      for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {


        G4VProcess const* proc = (*pVector)[j];
        G4String processName = proc->GetProcessName();

        //artificially attach generic name for parallel tracking worlds (only scoring planes live there)
        if (proc->GetProcessType() == G4ProcessType::fParallel) processName=G4String("Transportation");

        nUnknownProcesses+=insertIfNotFound(processName,particleName);


        if ( processName == gammaGP ){
          // special case for G4GammaGeneralProcess as of Geant4 11
          // if the process name is GammaGeneralProc then one has to
          // insert all its component processes:
          // phot, Rayl, compt, conv, GammaToMuPair, photonNuclear
          nUnknownProcesses+=insertIfNotFound(G4String("phot"),particleName);
          nUnknownProcesses+=insertIfNotFound(G4String("Rayl"),particleName);
          nUnknownProcesses+=insertIfNotFound(G4String("compt"),particleName);
          nUnknownProcesses+=insertIfNotFound(G4String("conv"),particleName);
          nUnknownProcesses+=insertIfNotFound(G4String("GammaToMuPair"),particleName);
          nUnknownProcesses+=insertIfNotFound(G4String("photonNuclear"),particleName);
        }

        // we will artificially attach "mu2eFieldPropagator" to all charged particles
        if( (particle->GetPDGCharge())!=0.0 ) {
          nUnknownProcesses+=insertIfNotFound(G4String("mu2eFieldPropagator"),particleName);
        }
        // we will artificially attach "NoProcess" to all particles
        // it was introduced in Geant4 v11.0 as a process name only and only for accounting purposes
        nUnknownProcesses+=insertIfNotFound(G4String("NoProcess"),particleName);

      } // end loop over processes
    }   // end loop over particle table

    // we will artificially attach "NotSpecified" process to Unspecified particle
    nUnknownProcesses+=insertIfNotFound(G4String("NotSpecified"),G4String("Unspecified"));

    if (nUnknownProcesses > 0 ){
      throw cet::exception("RANGE")
        << "There was one or more phyics processes that are not in the ProcessCode enum.\n"
        << "Number of processes: " << nUnknownProcesses
        << "\nPlease extend to the enum to add these processes and recompile.\n";
    }

    std::vector<ProcessCode> mu2eCodes = ProcessCode::mu2eCodes();
    for ( size_t i=0; i<mu2eCodes.size(); ++i){
      ProcessCode code = mu2eCodes[i];
      G4String name = code.name();
      _allProcesses.insert(std::make_pair(name,ProcInfo(code.name(),code) ));
    }
    // printAll(cout);

  } // PhysicsProcessInfo::beginRun

  void PhysicsProcessInfo::endRun(){

      //NEED TO PUT IN SOME PRINT SUMMARY FOR SEQUENTIAL MODE OR A SINGLE SUMMARRY IN MT
      //MODE.  RIGHT NOW WE GET ONE SUMMARY FOR EACH TRHEAD.
      //printSummary(cout);
    // printAll(cout);
  }

  // helper function to insert process and add particle name to ProcInfo; return non 0 if unknown
  int  PhysicsProcessInfo::insertIfNotFound( G4String const& processName,
                                             G4String const& particleName){
    int nUnknownProcesses(0);
    map_type::iterator jj = _allProcesses.find(processName);
    // If not already in the map, then create it.
    if ( jj == _allProcesses.end() ){
      _longestName = (processName.size() > _longestName) ? processName.size() : _longestName;

      ProcessCode code = ProcessCode::findByName(processName);
      if ( code.id() == ProcessCode::unknown ){
        ++nUnknownProcesses;
        cout << "Physics process named: " << processName
             << " is not known to the ProcessCode enum."
             << endl;
      }

      pair<map_type::iterator,bool> result =
        _allProcesses.insert( std::make_pair(processName,ProcInfo(processName,code) ));
      jj = result.first;
    }
    // add the particle to ProcInfo
    ProcInfo& ojj = jj->second;
    ojj.particleNames.push_back(particleName);
    return nUnknownProcesses;
  }

  // This can possibly be sped up considerably by checking frequently occuring names first.
  ProcessCode PhysicsProcessInfo::findAndCount( G4String const& name ){

    map_type::iterator i = _allProcesses.find(name);
    if ( i == _allProcesses.end() ){
      throw cet::exception("RANGE")
        << "Could not find physics process in PhysicsProcessInfo.  : "
        << name
        << "\n";
    }

    // Protect against overflowing the counters on very long jobs.
    if(i->second.count < std::numeric_limits<size_t>::max()) {
      ++i->second.count;
    }

    return i->second.code;
  }

  void PhysicsProcessInfo::printAll ( std::ostream& os) const{
    os << "Number of registered processes: "
       << _allProcesses.size() << " "
       << endl;

    int tcsum(0);
    for ( map_type::const_iterator i=_allProcesses.begin();
          i != _allProcesses.end(); ++i ){
      ProcInfo const& oi        = i->second;
      os << "Process: "
         << setw(_longestName) << oi.procName << " "
         << setw(4)  << oi.code.id()       << " "
         << setw(10) << oi.count;
      for ( size_t i=0; i<oi.particleNames.size(); ++i){
        os << " " << oi.particleNames[i];
      }
      os << endl;
      tcsum += oi.count;
    }

    os << "Total count: "
       << setw(_longestName) << " "
       << "  "
       << setw(10) << tcsum
       << endl;

  }


  void PhysicsProcessInfo::printSummary ( std::ostream& os) const{
    os << "Number of registered processes: "
       << _allProcesses.size() << " "
       << endl;

    int tcsum(0);
    for ( map_type::const_iterator i=_allProcesses.begin();
          i != _allProcesses.end(); ++i ){
      ProcInfo const& oi        = i->second;
      os << "Process: "
         << setw(_longestName) << oi.procName  << " "
         << setw(4)            << oi.code.id() << " "
         << setw(10)           << oi.count
         << endl;
      tcsum += oi.count;
    }

    os << "Total count: "
       << setw(_longestName) << " "
       << "  "
       << setw(10) << tcsum
       << endl;

  }

}  // end namespace mu2e
