//
// Types used in low-level tracker electronics 
// Dave Brown (LBNL) 25 Jun 2017
//
#ifndef TrackerConditions_Types_hh
#define TrackerConditions_Types_hh
#include <array>
#include <vector>

namespace mu2e {
  namespace TrkTypes {
    // 2 ends of the straw readout
    enum End{cal=0,hv, nends};
    // separately describe the 2 analog paths
    enum Path{thresh=0,adc,npaths};
    // typedefs used throught tracker 
    typedef unsigned long TDCValues[nends];
    typedef unsigned short TOTValues[nends];
    typedef std::array<double,nends> TDCTimes;
    typedef std::vector<unsigned short> ADCWaveform;
    typedef std::vector<double> ADCVoltages;
    typedef std::vector<double> ADCTimes;
  }
}
#endif
