//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALTrkClusterResiduals_hh
#define RecoDataProducts_ExtMonFNALTrkClusterResiduals_hh

namespace mu2e {

  //================================================================
  class ExtMonFNALTrkClusterResiduals {
  public:

    double dx() const { return dx_;}
    double xpull() const { return xpull_; }
    double dy() const { return dy_;}
    double ypull() const { return ypull_; }

    ExtMonFNALTrkClusterResiduals(double x, double xp, double y, double yp)
      : dx_(x), xpull_(xp), dy_(y), ypull_(yp)
    {}


    // ctr for genreflex persistency
    ExtMonFNALTrkClusterResiduals()
      : dx_(), xpull_(), dy_(), ypull_()
    {}

  private:
    double dx_;
    double xpull_;
    double dy_;
    double ypull_;
  };

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALTrkClusterResiduals_hh */
