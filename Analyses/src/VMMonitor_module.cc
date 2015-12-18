// Print current size of virtual memory usage for the currently executing process.
//
// Rob Kutschke, 2015

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

namespace {

  // A subset of the information from /proc/pid/status
  struct VMInfo {
    VMInfo():
      vmPeak(0), vmSize(0), vmHWM(0), vmRSS(0){
    }

    bool allOK() const {
      return haveVMPeak && haveVMSize && haveVMHWM && haveRSS;
    }

    void set ( std::string const& name, long val){
      if ( name == "VmPeak" ){
        vmPeak     = val;
        haveVMPeak = true;
      } else if ( name == "VmSize" ){
        vmSize     = val;
        haveVMSize = true;
      } else if ( name == "VmHWM"  ) {
        vmHWM      = val;
        haveVMHWM  = true;
      } else if ( name == "VmRSS"  ){
        vmRSS      = val;
        haveRSS    = true;
      }else{
        throw std::runtime_error( "VMMonitor: unrecognized name: " + name );
      }

    }

    long vmPeak;  // Peak virtual size
    long vmSize;  // Current virtual size
    long vmHWM;   // Max RSS ( high water mark )
    long vmRSS;   // Current RSS (Resident Set Size)

    bool haveVMPeak = false;
    bool haveVMSize = false;
    bool haveVMHWM  = false;
    bool haveRSS    = false;
  };

  // A helper class to populate and object of type VMInfo and serve information.
  class ProcStatus {
  public:
    ProcStatus();

    VMInfo const& info() const { return _proc; }

  private:
    long   _pageSize   = 0;
    size_t _maxVM      = 0;
    size_t _maxRSS     = 0;
    bool   _everyEvent = false;    // Info on every event or just on event with an increase.
    VMInfo _proc;
    int    _MB         = 1048576;  // Number of bytes in one MB.

  };

  // helper function
  long parseLine ( std::string const& line,
                   std::string const& keyExpected, 
                   std::string const& unitExpected ){
    std::istringstream is(line);
    std::string key, unit;
    long val;
    is >> key >> val >> unit;
    if ( unit!=unitExpected ){
      throw std::runtime_error( "VMMonitor: cannot parse: " + line );
    }
    return val;
  }

  // The c'tor does all of the work.
  // Extra info from /proc and populate the class member data.
  ProcStatus::ProcStatus(){
    std::ifstream  proc("/proc/self/status");
    std::array<std::string,4> wanted{ "VmPeak", "VmSize", "VmHWM", "VmRSS"};
    while ( proc ){
      std::string line;
      getline( proc, line);
      if ( !proc ) return;
      for ( auto const& name: wanted ){
        if ( line.find(name) != std::string::npos ){
          long val = parseLine( line , name, "kB");
          _proc.set(name,val);
          break;
        }
      }
    }
  }

} // end anonymous namespace

namespace mu2e {

  class VMMonitor : public art::EDAnalyzer {
  public:
    explicit VMMonitor(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
    long   _pageSize   = 0;        // Memory page size, in bytes.
    long   _kSize      = 1024;     // 1 kByte, in bytes.
    size_t _maxVM      = 0;        // Peak virtual size.
    size_t _maxRSS     = 0;        // Maximum resident set size.
    int    _fuzz       = 0;        // Only report large increases.
    bool   _everyEvent = false;    // Info on every event or just on event with an increase.
    int    _MB         = 1048576;  // Number of bytes in one MB.
  };

} // namespace mu2e

mu2e::VMMonitor::VMMonitor(const fhicl::ParameterSet& pset)
  : art::EDAnalyzer(pset),
    _fuzz(pset.get<int>("fuzz",0)),
    _everyEvent(pset.get<bool>("everyEvent",false))
{
  _pageSize = sysconf(_SC_PAGESIZE);
  if ( _pageSize == 0 ){
    throw std::runtime_error("Page size is zero!\n");
  }
}

void mu2e::VMMonitor::analyze(const art::Event& event) {

  // Extract virtual memory information from /proc/self/status
  // Convert to MBytes from kBytes
  ProcStatus vm;
  size_t vmPeak = (vm.info().vmPeak * _kSize)/_MB;
  size_t vmHWM  = (vm.info().vmHWM  * _kSize)/_MB;
  size_t vmSize = (vm.info().vmSize * _kSize)/_MB;
  size_t vmRSS  = (vm.info().vmRSS  * _kSize)/_MB;

  long deltaVM  = vmPeak - _maxVM;
  long deltaHWM = vmHWM  - _maxRSS;

  // Choose what to print, which may be nothing.
  int index = -1;
  if ( vmPeak > _maxVM+_fuzz ) {
    index = 0;
  } else if ( _everyEvent ){
    index = 1;
  }

  // Print it.
  static const std::array<std::string,2> blurb = {
    "VMMonitor: increasing VM or RSS at event: ",
    "VMMonitor: report VM and RSS at event:    "
  };
  if ( index > -1 && index < int(blurb.size()) ) {
    std::cout << blurb[index] << event.id()
              << " VM (MB): "   << vmPeak
              << " Delta(MB): " << deltaVM
              << " RSS(MB): "   << vmHWM
              << " Delta(MB): " << deltaHWM
              << " | Current VM(MB): " << vmSize
              << " Current RSS(MB): "  << vmRSS
              << std::endl;
  }

  // Increment peak counters
  _maxVM  = std::max(vmPeak,_maxVM);
  _maxRSS = std::max(vmHWM,_maxRSS);
}

void mu2e::VMMonitor::endJob() {
  std::cout << "VMMonitor: Peak virtual size: " << _maxVM  << std::endl;
  std::cout << "VMMonitor: Peak RSS:          " << _maxRSS << std::endl;
}

DEFINE_ART_MODULE(mu2e::VMMonitor);
