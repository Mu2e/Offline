///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_inc_Pzz_hh__
#define __CalPatRec_inc_Pzz_hh__

namespace CalPatRec {
  struct Pzz_t {
    int                              fID = 0;         // 3*face+panel, for pre-calculating overlaps
    double                           wx = 0.;         // direction cosines of the wire, all wires are assumed parallel
    double                           wy = 0.;
    double                           nx = 0.;         // direction cosines of the normal to the wires, pointing outwards
    double                           ny = 0.;
    float                            z = 0.f;         // Z-coordinate of the face
  };
}

#endif
