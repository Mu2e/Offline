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
    int  zface;                       // z-ordered face (so far, random, could've calculated, *TODO*)
    int  panel;                       // "geo" panel index
  };

  std::initializer_list<TrkPanelMap_t> TrkPanelMap_data = {
//-----------------------------------------------------------------------------
// station 0
// Z-ordering: plane 21: (101 | 219 | 235) outward  0
//                       (253 | 213 | 247) inward   1
//             plane 25: (248 | 276 | 224) inward   2
//                       (262 | 261 | 273) outward  3
//  mn_id  dtc lnk stn pln zf pnl
//-----------------------------------------------------------------------------
    { 261,  0,  0,  0,  0,  3, 0 },
    { 248,  0,  1,  0,  0,  2, 1 },
    { 224,  0,  2,  0,  0,  2, 2 },
    { 262,  0,  3,  0,  0,  3, 3 },
    { 273,  0,  4,  0,  0,  3, 4 },
    { 276,  0,  5,  0,  0,  2, 5 },
      
    { 253,  1,  0,  0,  1,  1, 0 },
    { 101,  1,  1,  0,  1,  0, 1 },
    { 219,  1,  2,  0,  1,  0, 2 },
    { 213,  1,  3,  0,  1,  1, 3 },
    { 235,  1,  4,  0,  1,  0, 4 },
    { 247,  1,  5,  0,  1,  1, 5 }
  };
};
#endif
