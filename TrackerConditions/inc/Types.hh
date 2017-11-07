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
//    typedef std::vector<unsigned short> ADCWaveform;
//    typedef std::array<unsigned short,16> ADCWaveform;
    typedef struct ADCWaveform {
      unsigned short& operator[] ( size_t index) { return _vals[index]; }
      unsigned short const& operator[] ( size_t index) const { return _vals[index]; }
      unsigned short _vals[NADC];
      constexpr size_t size() const { return NADC; }
      unsigned short& at( size_t index) { return _vals[index]; }
      unsigned short const& at( size_t index) const { return _vals[index]; }
      constexpr unsigned short* begin() { return &_vals[0]; }
      constexpr unsigned short* end() { return &_vals[NADC]; }
      constexpr const unsigned short* begin() const { return &_vals[0]; }
      constexpr const unsigned short* end() const { return &_vals[NADC]; }
    } ADCWaveform;
    typedef std::vector<float> ADCVoltages;
    typedef std::vector<float> ADCTimes;
  }
}
#endif
