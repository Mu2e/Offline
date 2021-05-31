// Print current size of virtual memory usage for the currently executing process.
//
// Rob Kutschke, 2015

#include "GeneralUtilities/inc/VMInfo.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class VMMonitor : public art::EDAnalyzer {
  public:
    explicit VMMonitor(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
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
{}

void mu2e::VMMonitor::analyze(const art::Event& event) {

  // Extract virtual memory information from /proc/self/status
  // Convert to MBytes from kBytes
  VMInfo vm;
  size_t vmPeak = (vm.vmPeak * _kSize)/_MB;
  size_t vmHWM  = (vm.vmHWM  * _kSize)/_MB;
  size_t vmSize = (vm.vmSize * _kSize)/_MB;
  size_t vmRSS  = (vm.vmRSS  * _kSize)/_MB;

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
