#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eStraightTrack.h"

using namespace mu2e;
namespace mu2e{

  TEveMu2eStraightTrack::TEveMu2eStraightTrack(){}

  /*------------Function to access position and direction:-------------*/
  XYZVectorF TEveMu2eStraightTrack::GetPositon(){ return fCosmicTrackSeed_->_track.MinuitEquation.Pos; }
  XYZVectorF TEveMu2eStraightTrack::GetDirection(){ return fCosmicTrackSeed_->_track.MinuitEquation.Dir; }

}
