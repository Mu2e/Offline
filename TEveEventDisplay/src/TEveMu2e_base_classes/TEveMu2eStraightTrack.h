#ifndef TEveMu2eStraightTrack_h
#define TEveMu2eStraightTrack_h

#include <TObject.h>
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include <TEveStraightLineSet.h>

namespace mu2e {
  class  TEveMu2eStraightTrack: public TEveStraightLineSet{

    public:
      #ifndef __CINT__
      explicit TEveMu2eStraightTrack();
      virtual ~TEveMu2eStraightTrack(){};

      CosmicTrackSeed* fCosmicTrackSeed_;

      void DrawStraightTrack();
      XYZVectorF GetPositon();
      XYZVectorF GetDirection();
      #endif
      ClassDef( TEveMu2eStraightTrack, 0);
  };
}
#endif
