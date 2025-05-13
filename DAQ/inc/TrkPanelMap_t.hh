///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __Offline_DAQ_inc_TrkPanelMap_t__
#define __Offline_DAQ_inc_TrkPanelMap_t__

namespace {
  struct TrkPanelMap_t {
    int  mnid;                        // 101='MN101' etc
    int  dtc;
    int  link;
    int  station;
    int  plane;                       // now: only one stations, plane = 0 or 1
    int  panel;                       // "geo" panel index
    int  zface;                       // z-ordered face (so far, random, could've calculated, *TODO*)
  };

  std::initializer_list<TrkPanelMap_t> TrkPanelMap_data = {
//-----------------------------------------------------------------------------
// station 0
// Z-ordering: plane 0 (21): (101 | 219 | 235) outward  0  
//                           (253 | 213 | 247) inward   1
//             plane 1 (25): (248 | 276 | 224) inward   2
//                           (262 | 261 | 273) outward  3
//  mn_id  dtc lnk stn pln pnl zf
//-----------------------------------------------------------------------------
    { 261, 18,  0,  0,  1,  3, 3},
    { 248, 18,  1,  0,  1,  4, 2},
    { 224, 18,  2,  0,  1,  0, 2},
    { 262, 18,  3,  0,  1,  5, 3},
    { 273, 18,  4,  0,  1,  1, 3},
    { 276, 18,  5,  0,  1,  2, 2},
                              
    { 253, 19,  0,  0,  0,  4, 1},
    { 101, 19,  1,  0,  0,  3, 0},
    { 219, 19,  2,  0,  0,  5, 0},
    { 213, 19,  3,  0,  0,  0, 1},
    { 235, 19,  4,  0,  0,  1, 0},
    { 247, 19,  5,  0,  0,  2, 1}
  };
};
#endif
