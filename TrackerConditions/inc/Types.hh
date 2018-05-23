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
    static constexpr size_t NADC = 16; // number of ADC samples including presamples
    static constexpr size_t NENDS = 2; // number of straw ends
    // typedefs used for tracker data structures 
    typedef std::array<uint16_t,NENDS> TDCValues;
    typedef std::array<uint16_t,NENDS> TOTValues;
    typedef std::array<float,NENDS> TDCTimes;
    typedef std::array<float,NENDS> TOTTimes;
    typedef std::array<uint16_t,NADC> ADCWaveform;
// transient types used in simulation
    typedef std::vector<float> ADCVoltages;
    typedef std::vector<float> ADCTimes;
  }
}
#endif
