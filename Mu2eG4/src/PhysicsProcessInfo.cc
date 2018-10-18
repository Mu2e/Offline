//
// Information about physics processes.
//
// $Id: PhysicsProcessInfo.cc,v 1.8 2012/12/13 23:57:58 genser Exp $
// $Author: genser $
// $Date: 2012/12/13 23:57:58 $
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
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

// G4 includes
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

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
        G4String const& name  = proc->GetProcessName();

        map_type::iterator jj = _allProcesses.find(name);

        // If not already in the map, then create it.
        if ( jj == _allProcesses.end() ){
          _longestName = (name.size() > _longestName) ? name.size() : _longestName;

          ProcessCode code = ProcessCode::findByName(name);
          if ( code.id() == ProcessCode::unknown ){
            ++nUnknownProcesses;
            cout << "Physics process named: " << name
                 << " is not known to the ProcessCode enum."
                 << endl;
          }

          pair<map_type::iterator,bool> result = _allProcesses.insert(
                 std::make_pair(name,ProcInfo(name,code) ));
          jj = result.first;

        }
        // we will artificially attach "FieldPropagator" to all particles
        // fixme : do it only for charged particles; factorize the code

        G4String const pname = G4String("FieldPropagator");
        jj = _allProcesses.find(pname);
        // If not already in the map, then create it.
        if ( jj == _allProcesses.end() ){
          _longestName = (pname.size() > _longestName) ? pname.size() : _longestName;

          ProcessCode code = ProcessCode::findByName(pname);
          if ( code.id() == ProcessCode::unknown ){
            ++nUnknownProcesses;
            cout << "Physics process named: " << pname
                 << " is not known to the ProcessCode enum."
                 << endl;
          }

          pair<map_type::iterator,bool> result = _allProcesses.insert(
                 std::make_pair(pname,ProcInfo(pname,code) ));
          jj = result.first;
        }

        // Add a particle to the
        ProcInfo& ojj = jj->second;
        ojj.particleNames.push_back(particleName);

      } // end loop over processes
    }   // end loop over particle table

    // we will artificially attach "NotSpecified" process to Unspecified particle
    // fixme factorize the code

    G4String pname = G4String("NotSpecified");
    _longestName = (pname.size() > _longestName) ? pname.size() : _longestName;

    ProcessCode pcode = ProcessCode::findByName(pname);
    if ( pcode.id() == ProcessCode::unknown ){
      ++nUnknownProcesses;
      cout << "Physics process named: " << pname
           << " is not known to the ProcessCode enum."
           << endl;
    }

    pair<map_type::iterator,bool> result =
      _allProcesses.insert(std::make_pair(pname,ProcInfo(pname,pcode)));
    map_type::iterator resultFirst = result.first;

    // Add Unspecified to the particleNames vector
    ProcInfo& pinfo = resultFirst ->second;
    pinfo.particleNames.push_back(G4String("Unspecified"));

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
