#ifndef Mu2eG4_PhysicsProcessInfo_hh
#define Mu2eG4_PhysicsProcessInfo_hh
//
// Information about physics processes.
//
//
// Original author Rob Kutschke
//

#include <iosfwd>
#include <set>
#include <map>
#include <string>

#include "MCDataProducts/inc/ProcessCode.hh"

#include "G4String.hh"

namespace mu2e {

  class PhysicsProcessInfo {

  public:

    PhysicsProcessInfo();

    // Accept compiler written d'tor, copy c'tor an assignment operator.

    // Build _allProcesses.  Valid until the end of the run.
    void beginRun();

    // Clear information at end of run and write summary printout.
    void endRun();

    // Locate a process by its name, return the corresponding process code and
    // increment the counter.
    ProcessCode findAndCount( G4String const& name );
      
    void printAll ( std::ostream& os) const;
    void printSummary ( std::ostream& os) const;

    // Information about one physics process.
    struct ProcInfo{
      ProcInfo():procName(""),particleNames(),code(),count(0){}
      ProcInfo( std::string aprocName,
             ProcessCode acode ):
        procName(aprocName),
        particleNames(),
        code(acode),
        count(0){}

      std::string procName;
      std::vector<std::string> particleNames;
      ProcessCode code;
      size_t count;
    };

  private:

    // The information about all of the physics processes.
    typedef std::map<G4String,ProcInfo> map_type;
    map_type _allProcesses;

    // The length of the longest name; for formatting printed output.
    size_t _longestName;

  };

} // end namespace mu2e

#endif /* Mu2eG4_PhysicsProcessInfo_hh */

