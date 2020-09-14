#ifndef TEveMu2eStraightTrack_h
#define TEveMu2eStraightTrack_h

#include <TObject.h>
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "DataProducts/inc/XYZVec.hh"
#include <TEveStraightLineSet.h>

namespace mu2e {
  class  TEveMu2eStraightTrack: public TEveStraightLineSet{
      CosmicTrackSeed* fCosmicTrackSeed;
    public:
      #ifndef __CINT__
      explicit TEveMu2eStraightTrack(){};
      virtual ~TEveMu2eStraightTrack(){};
      #endif
      void DrawStraightTrack();
      XYZVec GetPositon(){ return fCosmicTrackSeed->_track.MinuitEquation.Pos; } ;
      XYZVec GetDirection(){ return fCosmicTrackSeed->_track.MinuitEquation.Dir; };
      ClassDef( TEveMu2eStraightTrack, 0);
  };
}
#endif
