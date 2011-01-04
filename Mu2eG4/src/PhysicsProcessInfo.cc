//
// Information about physics processes.
//
// $Id: PhysicsProcessInfo.cc,v 1.2 2011/01/04 22:10:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/01/04 22:10:23 $
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
#include "FWCore/Utilities/interface/Exception.h"

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
      G4ProcessVector const* pVector      = pmanager->GetProcessList();
      G4String               particleName = particle->GetParticleName();

      // Loop over all processes.
      for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) {

        G4VProcess const* proc = (*pVector)[j];
        G4String const& name  = proc->GetProcessName();

        map_type::iterator jj = _allProcesses.find(name);
        
        // If not already in the map, then create it.
        if ( jj == _allProcesses.end() ){
          _longestName = (name.size() > _longestName) ? name.size() : _longestName;

          ProcessCode code = ProcessCode::findByName(name);
          if ( code.id() == ProcessCode::unknown ){
            ++nUnknownProcesses;
            cout << "Phyics process named: " << name
                 << " is not known to the ProcessCode enum."
                 << endl;
          }

          pair<map_type::iterator,bool> result = _allProcesses.insert( 
                 std::make_pair(proc->GetProcessName(),ProcInfo(name,code) ));
          jj = result.first;

        }

        // Add a particle to the 
        ProcInfo& ojj = jj->second;
        ojj.particleNames.push_back(particleName);

      } // end loop over processes
    }   // end loop over particle table

    if (nUnknownProcesses > 0 ){
      throw cms::Exception("RANGE")
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

  } // PhysicsProcessInfo::beginRun 

  void PhysicsProcessInfo::endRun(){
    printSummary(cout);
  }

  // This can possibly be sped up considerably by checking frequently occuring names first.
  ProcessCode PhysicsProcessInfo::findAndCount( G4String const& name ){
    map_type::iterator i = _allProcesses.find(name);
    if ( i == _allProcesses.end() ){
      throw cms::Exception("RANGE")
        << "Could not find physics process in PhysicsProcessInfo.  : "
        << name
        << "\n";
    }

    // Protect against overflowing the counters on very long jobs.
    static const int bigNumber(1000000000);
    static int counter(0);
    if ( counter < bigNumber ){
      ++counter;
      ++(i->second.count);
    }

    return i->second.code;
  }

  void PhysicsProcessInfo::printAll ( std::ostream& os) const{
    os << "Number of registered processes: "
       << _allProcesses.size() << " "
       << endl;

    for ( map_type::const_iterator i=_allProcesses.begin();
          i != _allProcesses.end(); ++i ){
      ProcInfo const& oi        = i->second;
      os << "Process: "
         << setw(_longestName) << oi.procName << " "
         << oi.code            << " "
         << oi.count;
      for ( size_t i=0; i<oi.particleNames.size(); ++i){
        os << " " << oi.particleNames[i];
      }
      os << endl;
    }
  }


  void PhysicsProcessInfo::printSummary ( std::ostream& os) const{
    os << "Number of registered processes: "
       << _allProcesses.size() << " "
       << endl;
    
    for ( map_type::const_iterator i=_allProcesses.begin();
          i != _allProcesses.end(); ++i ){
      ProcInfo const& oi        = i->second;
      os << "Process: "
         << setw(_longestName) << oi.procName  << " "
         << setw(4)            << oi.code.id() << " "
         << oi.count
         << endl;
      
    }
  }

  

}  // end namespace mu2e
