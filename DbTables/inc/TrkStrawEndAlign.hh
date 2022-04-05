#ifndef DbTables_TrkStrawEndAlign_hh
#define DbTables_TrkStrawEndAlign_hh
//
// Struct describing displacements of the straw and wire ends
//
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  struct TrkStrawEndAlign {
    using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO
    TrkStrawEndAlign(int index, StrawId const& id,
        float wire_cal_dV, float wire_cal_dW,
        float wire_hv_dV, float wire_hv_dW,
        float straw_cal_dV, float straw_cal_dW,
        float straw_hv_dV, float straw_hv_dW) :
      _wire_cal_dV(wire_cal_dV), _wire_cal_dW(wire_cal_dW),
      _wire_hv_dV(wire_hv_dV), _wire_hv_dW(wire_hv_dW),
      _straw_cal_dV(straw_cal_dV), _straw_cal_dW(straw_cal_dW),
      _straw_hv_dV(straw_hv_dV), _straw_hv_dW(straw_hv_dW) {}
    StrawId const& id() const { return _id; }
    xyzVec wireDeltaUVW(StrawEnd::End end) const { return end == StrawEnd::cal ? xyzVec(0.0,_wire_cal_dV,_wire_cal_dW) : xyzVec(0.0,_wire_hv_dV,_wire_hv_dW); }
    xyzVec wireDeltaUVW(StrawEnd const& end) const { return wireDeltaUVW(end.end()); }
    xyzVec strawDeltaUVW(StrawEnd::End end) const { return end == StrawEnd::cal ? xyzVec(0.0,_straw_cal_dV,_straw_cal_dW) : xyzVec(0.0,_straw_hv_dV,_straw_hv_dW); }
    xyzVec strawDeltaUVW(StrawEnd const& end) const { return strawDeltaUVW(end.end()); }

    int _index;
    StrawId _id;
    float _wire_cal_dV, _wire_cal_dW; // wire cal end displacements from nominal
    float _wire_hv_dV, _wire_hv_dW; // wire hv end displacements from nominal
    float _straw_cal_dV, _straw_cal_dW; // straw cal end displacements from nominal
    float _straw_hv_dV, _straw_hv_dW; // straw hv end displacements from nominal
  };


}
#endif
