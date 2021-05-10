#ifndef RecoDataProducts_LumiInfo_hh
#define RecoDataProducts_LumiInfo_hh
#include <float.h>
#include <math.h>
//
// This object is a mock-up to contain the variables to monitor the luminosity
// The exact structure is still under discussion. (now just an unsigned short [5])
//
// To save space we might aim to have a unit16_t * lumi = new unit16_t; with
// lumi[0] : 8 bit for the version, 8 bit for the number of entries
// lumi[i] : 16 bit to indicate the i-number for the lumi (proton TC, calo energy etc)
//
// bvitali May 2021


namespace mu2e {

  struct LumiInfo {
    
    LumiInfo(){
      _lumi[0] = 0;
      _lumi[1] = 0;
      _lumi[2] = 0;
      _lumi[3] = 0;
      _lumi[4] = 0;
    }    

    unsigned short _lumi[5];

  };
}

#endif


