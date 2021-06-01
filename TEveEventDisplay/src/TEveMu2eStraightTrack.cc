#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eStraightTrack.h"

using namespace mu2e;
namespace mu2e{

  TEveMu2eStraightTrack::TEveMu2eStraightTrack(){};
  
  /*------------Function to access position and direction:-------------*/
  XYZVec TEveMu2eStraightTrack::GetPositon(){ return fCosmicTrackSeed_->_track.MinuitEquation.Pos; } ;
  XYZVec TEveMu2eStraightTrack::GetDirection(){ return fCosmicTrackSeed_->_track.MinuitEquation.Dir; };
  
}
