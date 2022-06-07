//
// Helper class to hold functions formerly in RecoDataProducts/inc/CosmicTrack.hh but which
// use classes from Mu2eUtilities.
//
#include "Offline/Mu2eUtilities/inc/CosmicTrackUtils.hh"

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

namespace mu2e {

  std::tuple <double, double, double, double, double, double> KinKalTrackParams( CosmicTrack const& trk){
    XYZVectorF zpos(0.,0.,0);
    XYZVectorF zdir(0.,0.,1.);
    XYZVectorF pos0(trk.MinuitParams.A0, 0, trk.MinuitParams.B0);
    XYZVectorF dir(trk.MinuitParams.A1, -1, trk.MinuitParams.B1);

    std::tuple <double,double, double, double, double, double> info;
    TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(pos0, dir, zpos, zdir);
    XYZVectorF POCA = PCA.point1()-PCA.point2();
    double DOCA = PCA.dca();
    double amsign = copysign(1.0, -(zdir.Cross(POCA)).Dot(dir));

    double d0 = amsign*DOCA;
    double phi0 = dir.Phi();
    double z0 = PCA.point1().Z();
    double cost = dir.Z();
    double t0 = trk.MinuitParams.T0; //TODO
    double mom = 1.0;//TODO
    info = std::make_tuple(d0,phi0,z0,cost, t0, mom);
    return info;
  }

}
