#ifndef __CalPatRec_inc_ChannelID_hh__
#define __CalPatRec_inc_ChannelID_hh__

namespace CalPatRec {
  struct ChannelID {
    int Station;
    int Plane;
    int Face;
    int Panel;

    static void orderID  (ChannelID* X, ChannelID* Ordered);
    static void deOrderID(ChannelID* X, ChannelID* Ordered);

    static void testOrderID  ();
    static void testdeOrderID();
  };
}

#endif
