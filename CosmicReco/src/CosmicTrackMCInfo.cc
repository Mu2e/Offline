#include "CosmicReco/inc/CosmicTrackMCInfo.hh"

CosmicTrackMCInfo::CosmicTrackMCInfo(){
     
     TrueFitEquation.Pos.SetXYZ(0,0,0);
     TrueFitEquation.Dir.SetXYZ(0,0,0);
     TrueTrackCoordSystem._XDoublePrime.SetXYZ(0,0,0);
     TrueTrackCoordSystem._YDoublePrime.SetXYZ(0,0,0);
     TrueTrackCoordSystem._ZPrime.SetXYZ(0,0,0);
     RawTrueParams.A0 = 0.;
     RawTrueParams.A1 = 0.;
     RawTrueParams.B0 =0;
     RawTrueParams.B1 = 0.;
	
} 


