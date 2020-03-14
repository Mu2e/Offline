#include "CosmicReco/inc/CosmicTrackMCInfo.hh"

CosmicTrackMCInfo::CosmicTrackMCInfo(){
     TrueTheta =0;
     TruePhi =0;
     
     TruePhiSIM = 0;
     TrueThetaSIM = 0;
     TrueMomentum = 0;

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


