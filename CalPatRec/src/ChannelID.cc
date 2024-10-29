///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include "Offline/CalPatRec/inc/ChannelID.hh"

namespace CalPatRec {
//-----------------------------------------------------------------------------
    void ChannelID::orderID(ChannelID* X, ChannelID* O) {
      if (X->Panel % 2 == 0) X->Face = 0;
      else                   X->Face = 1; // define original face

      O->Station = X->Station;            // stations already ordered
      O->Plane   = X->Plane;              // planes already ordered. Not using them

      if (X->Station % 2 == 0) {
//-----------------------------------------------------------------------------
// even stations : 0,2,4...
//-----------------------------------------------------------------------------
        if (X->Plane == 0) O->Face = 1 - X->Face;
        else               O->Face = X->Face + 2;
      }
      else {
//-----------------------------------------------------------------------------
// odd stations : 1,3,5...
//-----------------------------------------------------------------------------
        if (X->Plane == 0) O->Face = X->Face;
        else               O->Face = 3 - X->Face; // order face
      }

      O->Panel = int(X->Panel/2);                // 3 panels per face, face0: 024
    }

//-----------------------------------------------------------------------------
  void ChannelID::deOrderID(ChannelID* X, ChannelID* O) {

    X->Station = O->Station;
    X->Plane   = O->Plane;

    if(O->Station % 2 ==  0) {
      if(O->Plane == 0) X->Face = 1 - O->Face;
      else X->Face = O->Face - 2;
    }
    else {
      if(O->Plane == 0) X->Face = O->Face;
      else X->Face = 3 - O->Face;
    }

    if(X->Face == 0) X->Panel = O->Panel * 2;
    else X->Panel = 1 + (O->Panel * 2);
  }

//-----------------------------------------------------------------------------
// testOrderID & testdeOrderID not used in module, they are only used to make sure
// OrderID and deOrderID work as intended
//-----------------------------------------------------------------------------
  void ChannelID::testdeOrderID() {

    ChannelID x, o;

    for (int s=0; s<2; ++s) {
      for (int f=0; f<4; ++f) {
        for (int pa=0; pa<3; ++pa) {
          o.Station          = s;
          o.Face             = f;
          if (f < 2) o.Plane = 0;
          else       o.Plane = 1;
          o.Panel            = pa;

          deOrderID(&x, &o);

          printf(" testdeOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i)",
                 x.Station,x.Plane,x.Face,x.Panel);
          printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i)\n",
                 o.Station,o.Plane,o.Face,o.Panel);
        }
      }
    }
  }

//-----------------------------------------------------------------------------
  void ChannelID::testOrderID() {

    ChannelID x, o;

    for (int s=0; s<2; ++s) {
      for (int pl=0; pl<2; ++pl) {
        for (int pa=0; pa<6; ++pa) {
          x.Station = s;
          x.Plane   = pl;
          x.Panel   = pa;
          orderID(&x, &o);
          printf(" testOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i)",
                 x.Station, x.Plane, x.Face, x.Panel);
          printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i)\n",
                 o.Station, o.Plane, o.Face, o.Panel);
        }
      }
    }
  }

}
