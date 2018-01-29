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
    static constexpr size_t NADC = 16; // this is also a parameter in StrawElectronics, FIXME!
    // separately describe the 2 analog paths
    enum Path{thresh=0,adc,npaths};
    // typedefs used for tracker data structures 
    typedef std::array<uint16_t,nends> TDCValues;
    typedef std::array<uint16_t,nends> TOTValues;
    typedef std::array<float,nends> TDCTimes;
    typedef std::array<float,nends> TOTTimes;
    typedef std::array<uint16_t,NADC> ADCWaveform;
// transient types used in simulation
    typedef std::vector<float> ADCVoltages;
    typedef std::vector<float> ADCTimes;
  }
}
#endif
