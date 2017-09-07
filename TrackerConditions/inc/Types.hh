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
    // typedefs used for tracker data structures 
//    typedef std::array<unsigned long,nends> TDCValues;
//    typedef std::array<unsigned short,nends> TOTValues;
//    typedef std::array<float,nends> TDCTimes;
    typedef struct TDCValues { 
      unsigned long& operator[] ( size_t index) { return _vals[index]; }
      unsigned long const& operator[] ( size_t index) const { return _vals[index]; }
      unsigned long _vals[nends];
    } TDCValues;
    typedef struct TOTValues { 
      unsigned short& operator[] ( size_t index) { return _vals[index]; }
      unsigned short const& operator[] ( size_t index) const { return _vals[index]; }
      unsigned short _vals[nends];
    } TOTValues;
    typedef struct TDCTimes { 
      float& operator[] ( size_t index) { return _vals[index]; }
      float const& operator[] ( size_t index) const { return _vals[index]; }
      float _vals[nends];
    } TDCTimes;
    typedef struct TOTTimes { 
      float& operator[] ( size_t index) { return _vals[index]; }
      float const& operator[] ( size_t index) const { return _vals[index]; }
      float _vals[nends];
    } TOTTimes;
    typedef std::vector<unsigned short> ADCWaveform;
    typedef std::vector<float> ADCVoltages;
    typedef std::vector<float> ADCTimes;
  }
}
#endif
