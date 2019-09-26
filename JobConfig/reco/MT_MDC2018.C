//
// Generate random misalignemnts for MDC2018.  Execute as:
//  root -l JobConfig/reco/MT_MDC2018.C
//
#include "JobConfig/reco/MisalignTracker.C+"
#include "math.h"
{
  double twist = 5; // 5 mm at the outer radius over the length of the tracker (~1 mrad)
  double skew = 5; // 5mm displacement over the length of the tracker
  double squeeze = 5; // 5mm over the length of the tracker
  double anglesig = 0.0007; // 0.7 mrad, ~500 um at the outer edge
  double possig = 0.5; // 500 um precision
  MisalignTracker(twist,skew,squeeze,anglesig,possig,"MT_MDC2018.txt")
}
