#ifndef GeneralUtilities_MVAStatus_hh
#define GeneralUtilities_MVAStatus_hh
#include <iostream>
namespace mu2e {
  struct MVAStatus {
    enum mvastat {unset=0,filled,calculated,failed};
    int16_t _status;
    MVAStatus() : _status(unset){}
    MVAStatus(mvastat status) : _status(status) {}
    operator int16_t() { return _status; }
  };
}
#endif
