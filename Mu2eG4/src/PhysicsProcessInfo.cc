//
// Information about physics processes.
//
// $Id: PhysicsProcessInfo.cc,v 1.1 2010/12/17 22:05:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/17 22:05:56 $
//
// Original author Rob Kutschke
//
// Notes
// 1) About the iterator class:
//      G4ParticleTable::G4PTblDicIterator
//    The G4ParticleTable class internally contains a 
//      std::map<G4String,ParticleDefinition*>
//    Users of the class need to be able to loop over the map.  To provide 
//    this functionality the designers of the class chose to provide a 
//    custom written iterator.  ( Current best practice would be to provide
//    a typedef to the underlying stl iterator type.  It is possible that
//    this design was frozen before std::map was robust and that this class
//    provides a work around that once had been necessary. )
//
//    The rules for using this iterator are:
//     1) You must call reset() before starting the loop.
//        If you forget you will usually, but not always, start the iteration
//        past the end of the map, which produces undefined behaviour.
//     2) The call to operator ():
//          a) on the first call after reset, it will initialize an internal 
//             iterator to point at the first element of the map.
//          b) on subsequent calls it will increment the internal iterator.
//          c) On all calls it will return true(false) if the interal iterator
//             is !=(==) to the .end() of the map.
//     3) At any time, a call to key() or to value() will return the information
//        pointed to by the interally held iterator.
//
//     A particularly bizarre feature is that the G4ParticleTable class holds
//     its own iterator.  When you ask for an iterator, you actually get a pointer
//     to this one iterator.  Therefore nested loops over the map have totally
//     unexpected behaviour.

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

          StoppingCode code = StoppingCode::findByName(name);
          if ( code.id() == StoppingCode::unknown ){
            throw cms::Exception("RANGE")
              << "Phyics process named: " << name
              << " is not known to the StoppingCode enum.\n"
              << "Please extend to the enum to add this process.\n";
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

    std::vector<StoppingCode> mu2eCodes = StoppingCode::mu2eCodes();
    for ( size_t i=0; i<mu2eCodes.size(); ++i){
      StoppingCode code = mu2eCodes[i];
      G4String name = code.name();
      _allProcesses.insert(std::make_pair(name,ProcInfo(code.name(),code) ));
    }


  } // PhysicsProcessInfo::beginRun 

  void PhysicsProcessInfo::endRun(){
    printSummary(cout);
  }

  // This can possibly be sped up considerably by checking frequently occuring names first.
  StoppingCode PhysicsProcessInfo::findAndCount( G4String const& name ){
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
