//
// Types used in low-level tracker electronics
// Dave Brown (LBNL) 25 Jun 2017
//
#ifndef TrackerConditions_Types_hh
#define TrackerConditions_Types_hh
#include <array>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace mu2e {
  namespace TrkTypes {
    // 2 ends of the straw readout
    static constexpr size_t NADC_MIN = 3; // number of ADC samples included in main packet
    static constexpr size_t NADC_PERPACKET = 12; // number of ADC samples in each additional packet
    static constexpr size_t NENDS = 2; // number of straw ends
    // typedefs used for tracker data structures
    typedef uint32_t TDCValue;
    typedef uint16_t TOTValue;
    typedef uint16_t ADCValue;
    typedef std::array<TDCValue,NENDS> TDCValues;
    typedef std::array<TOTValue,NENDS> TOTValues;
    typedef std::array<float,NENDS> TDCTimes;
    typedef std::array<float,NENDS> TOTTimes;
    typedef std::vector<ADCValue> ADCWaveform;
// transient types used in simulation
    typedef std::vector<float> ADCVoltages;
    typedef std::vector<float> ADCTimes;
  }
}
#endif
