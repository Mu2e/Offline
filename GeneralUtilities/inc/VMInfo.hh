#ifndef Mu2eUtilities_VMInfo_hh
#define Mu2eUtilities_VMInfo_hh

// Information about virtual memory usage.
// A subset of the information from /proc/pid/status.
//
// Original Author Rob Kutschke, December 2020

namespace mu2e {

  struct VMInfo {
    VMInfo();

    long vmPeak = 0;  // Peak virtual size, in KiB
    long vmSize = 0;  // Current virtual size, in KiB
    long vmHWM  = 0;  // Max RSS ( Virtual Memory High Water Mark ), in KiB
    long vmRSS  = 0;  // Current Resident Set Size, in KiB

  }; // struct VMInfo

} // namespace mu2e

#endif // Mu2eUtilities_VMInfo_hh
