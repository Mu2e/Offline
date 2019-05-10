//
//  Simple object to keep track of original object counts before reconstruction compression
//  Original author: Dave Brown (LBNL) Feb 2019
//
#ifndef RecoDataProducts_RecoCount_HH
#define RecoDataProducts_RecoCount_HH
#include "DataProducts/inc/AHist.hh"

namespace mu2e {
 struct RecoCount {
    uint16_t _nstrawdigi; // original number of straw digis in the event   
    uint16_t _ncalodigi; // original number of calo digis in the event   
    uint16_t _ncrvdigi; // original number of crv digis in the event   
    uint16_t _nshfesel, _nshfrsel, _nshftsel, _nshfbkg, _nshftpk; // StrawHitFlag counts
    // histogram of straw hit times
    static constexpr size_t _nshtbins =70; // number of bins in the Straw Hit time histogram
    AHist<uint16_t, _nshtbins> _shthist;
    RecoCount() : _nstrawdigi(0), _ncalodigi(0),_ncrvdigi(0), 
    _nshfesel(0), _nshfrsel(0), _nshftsel(0), _nshfbkg(0), _nshftpk(0),
    _shthist(-25.0,1725.0) {}
  };
}
#endif
