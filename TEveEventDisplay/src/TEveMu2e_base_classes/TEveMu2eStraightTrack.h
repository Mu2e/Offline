#ifndef TEveMu2eStraightTrack_h
#define TEveMu2eStraightTrack_h

#include <TObject.h>
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wtype-limits"
//#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
//#pragma GCC diagnostic pop

#include <TEveStraightLineSet.h>

namespace mu2e {
  class  TEveMu2eStraightTrack: public TEveStraightLineSet{
      //CosmicTrackSeed* fCosmicTrackSeed;
    public:
      #ifndef __CINT__
      explicit TEveMu2eStraightTrack(){};
      virtual ~TEveMu2eStraightTrack(){};
      #endif
      void DrawStraightTrack();
      void GetPositon();
      void GetDirection();
      ClassDef( TEveMu2eStraightTrack, 0);
  };
}
#endif
