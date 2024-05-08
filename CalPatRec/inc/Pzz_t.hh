///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_inc_Pzz_hh__
#define __CalPatRec_inc_Pzz_hh__

namespace CalPatRec {
  struct Pzz_t {
    int                              fID;         // 3*face+panel, for pre-calculating overlaps
    double                           wx;          // direction cosines of the wire, all wires are assumed parallel
    double                           wy;
    double                           nx;          // direction cosines of the normal to the wires, pointing outwards
    double                           ny;
    float                            z;           // Z-coordinate of the face
  };
}

#endif
