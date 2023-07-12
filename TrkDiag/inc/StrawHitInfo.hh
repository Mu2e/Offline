//
// struct for hit diagnostics
#ifndef StrawHitInfo_hh
#define StrawHitInfo_hh
#include "Offline/DataProducts/inc/GenVector.hh"

namespace mu2e {
  struct StrawHitInfo {
    XYZVectorF _pos;
    float _edep = 0.0;
    float _time = 0.0;
    float _wdist = 0.0;
    float _wres = 0.0;
    float _tres = 0.0;
    float _chisq = 0.0;
    float _stdt = 0.0;
    float _dist = 0.0;
    int _plane = 0.0;
    int _panel = 0.0;
    int _layer = 0.0;
    int _straw = 0.0;
    bool _bkg = false;
    bool _bkgc = false;
    bool _isolated = false;
    bool _tsel = false;
    bool _esel = false;
    bool _rsel = false;
    bool _stereo = false;
    bool _tdiv = false;
    bool _strawxtalk = false;
    bool _elecxtalk = false;
    XYZVectorF _mcpos;
    int _mcpdg = 0;
    int _mcproc = 0;
    float _mcedep = 0.0;
    float _mctime = 0.0;
    float _mcht = 0.0;
    float _mcmom = 0.0;
    float _mctd = 0.0;
    int _prel = -1;
    int _mrel = -1;
  };
}
#endif
