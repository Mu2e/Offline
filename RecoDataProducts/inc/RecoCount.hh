//
//  Simple object to keep track of original object counts before reconstruction compression
//  Original author: Dave Brown (LBNL) Feb 2019
//
#ifndef RecoDataProducts_RecoCount_HH
#define RecoDataProducts_RecoCount_HH
namespace mu2e {
  struct RecoCount {
    uint16_t _nstrawdigi; // original number of straw digis in the event   
    uint16_t _ncalodigi; // original number of calo digis in the event   
    uint16_t _ncrvdigi; // original number of crv digis in the event   
    uint16_t _ncaloclust; // original number of calo clusters in the event
    uint16_t _nshfesel, _nshfrsel, _nshftsel, _nshfbkg, _nshftpk; // StrawHitFlag counts
    RecoCount() : _nstrawdigi(0), _ncalodigi(0),_ncrvdigi(0), _ncaloclust(0),
    _nshfesel(0), _nshfrsel(0), _nshftsel(0), _nshfbkg(0), _nshftpk(0) {}
  };
}
#endif
