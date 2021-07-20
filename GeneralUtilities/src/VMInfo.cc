#include "Offline/GeneralUtilities/inc/VMInfo.hh"

#include <array>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

namespace {

  // Helper function to check the units and return the value.
  long parseLine ( std::string const& line,
                   std::string const& unitExpected ){
    std::istringstream is(line);
    std::string key, unit;
    long val;
    is >> key >> val >> unit;
    if ( unit!=unitExpected ){
      throw std::runtime_error( "ProcStatus: cannot parse: " + line );
    }
    return val;
  }

} // end anonymous namespace

// The c'tor does all of the work.
mu2e::VMInfo::VMInfo(){

  // The information we want.
  std::array<std::string,4> wanted{ "VmPeak", "VmSize", "VmHWM", "VmRSS"};

  // Parse the information to get what we need.
  std::ifstream  proc("/proc/self/status");
  std::map<std::string,long> values;
  while ( proc ){
    std::string line;
    getline( proc, line);
    if ( !proc ) break;
    for ( auto const& name: wanted ){
      if ( line.find(name) != std::string::npos ){
        long val = parseLine( line,  "kB");
        values[name] = val;
        break;
      }
    }
  }

  vmPeak = values["VmPeak"];
  vmSize = values["VmSize"];
  vmHWM  = values["VmHWM"];
  vmRSS  = values["VmRSS"];
}
