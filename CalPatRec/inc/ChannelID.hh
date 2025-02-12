#ifndef __CalPatRec_inc_ChannelID_hh__
#define __CalPatRec_inc_ChannelID_hh__

namespace CalPatRec {
  struct ChannelID {
    int Station = 0;
    int Plane = 0;
    int Face = 0;
    int Panel = 0;

    static void orderID  (ChannelID* X, ChannelID* Ordered);
    static void deOrderID(ChannelID* X, ChannelID* Ordered);

    static void testOrderID  ();
    static void testdeOrderID();
  };
}

#endif
