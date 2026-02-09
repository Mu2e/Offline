// TrackSummary is a persistable class to communicate track
// reconstruction results to physics analyses.
//
// Andrei Gaponenko, 2014

#ifndef RecoDataProducts_TrackSummary_hh
#define RecoDataProducts_TrackSummary_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <ostream>

class TrkSimpTraj; // BaBar class

namespace mu2e {


  class TrackSummary {
  public:
    
    //================================================================
    class HelixParams {
    public:
      // See docdb-781 for parameter definitions
      double d0() const { return d0_; }
      double phi0() const { return phi0_; }
      double omega() const { return omega_; }
      double z0() const { return z0_; }
      double tanDip() const { return tanDip_; }

      // Picked from BTrk/BaBar/BTrk/TrkBase/include/HelixParams.hh
      enum ParIndex {d0Index=0, phi0Index, omegaIndex, z0Index, tanDipIndex, NHLXPRM};
      const CLHEP::HepSymMatrix& covariance() const { return covariance_; }

      //Not yet implemented - commented out until it is.
      //double parErr(ParIndex i);

      // Some derived quantities
      double dOut() const; // max distance to Z axis, opposite to d0()
      double radius() const; // of the helix
      double wavelength() const; // of the helix

      explicit HelixParams(const TrkSimpTraj& ltraj);

      // Default constructor is required by ROOT persistency
      HelixParams() : d0_(), phi0_(), omega_(), z0_(), tanDip_(), covariance_(NHLXPRM) {}

    private:
      double d0_;
      double phi0_;
      double omega_;
      double z0_;
      double tanDip_;
      CLHEP::HepSymMatrix covariance_;
    };
    
    //================================================================
    class TrackStateAtPoint {
    public:
      const HelixParams& helix() const { return helix_; }
      const CLHEP::Hep3Vector& position() const { return position_; }
      const CLHEP::Hep3Vector& momentum() const { return momentum_; }
      const CLHEP::HepSymMatrix& momentumCovariance() const { return momentumCovariance_; }

      double arrivalTime() const { return arrivalTime_; }
      double flightLength() const { return flightLength_; }

      // Derived quantities
      double momentumError() const;
      double costh() const; // momentum vector cosine to the Z axis

      TrackStateAtPoint(const HelixParams& h,
                        const CLHEP::Hep3Vector& pos,

                        const CLHEP::Hep3Vector& mom, const CLHEP::HepSymMatrix& momCov,
                        double arrivalTime, double flightLength)
        : helix_(h)
        , position_(pos), momentum_(mom)
        , momentumCovariance_(momCov)
        , arrivalTime_(arrivalTime), flightLength_(flightLength)
      {}

      // Default constructor is required by ROOT persistency
      TrackStateAtPoint() : arrivalTime_(), flightLength_() {}

    private:
      // Converting (pos,mom) to a helix depends on B field and we want
      // to decouple from that.  Therefore we store redundant info here.
      HelixParams helix_;
      CLHEP::Hep3Vector position_;
      CLHEP::Hep3Vector momentum_;
      CLHEP::HepSymMatrix momentumCovariance_;
      double arrivalTime_;
      double flightLength_;
    };

    //================================================================
    int fitstatus() const { return fitstatus_; }
    int charge() const { return charge_; }
    int nactive() const { return nactive_; }

    int ndof() const { return ndof_; }
    double chi2() const { return chi2_; }
    double fitcon() const;

    double t0() const { return t0_; }
    double t0Err() const { return t0Err_; }
    // flight length associated with t0.
    double flt0() const { return flt0_; }

    const std::vector<TrackStateAtPoint>& states() const { return states_; }

    TrackSummary(int fitstatus, int charge, int nactive,
                 int ndof, double chi2,
                 double t0, double t0Err, double flt0)
      : fitstatus_(fitstatus), charge_(charge), nactive_(nactive)
      , ndof_(ndof), chi2_(chi2)
      , t0_(t0), t0Err_(t0Err), flt0_(flt0)
    {}

    void addState(const TrackStateAtPoint& st);

    // Default constructor is required by ROOT persistency
    TrackSummary() : fitstatus_(), charge_(), nactive_(), ndof_(), chi2_(), t0_(), t0Err_(), flt0_() {}

  private:
    std::vector<TrackStateAtPoint> states_;

    int fitstatus_;
    int charge_;
    int nactive_;

    int ndof_;
    double chi2_;

    double t0_;
    double t0Err_;
    double flt0_;

  };

  //================================================================
  typedef  std::vector<TrackSummary>  TrackSummaryCollection;

  std::ostream& operator<<(std::ostream& os, const TrackSummary::TrackStateAtPoint& st);
  std::ostream& operator<<(std::ostream& os, const TrackSummary& sum);
  std::ostream& operator<<(std::ostream& os, const TrackSummaryCollection& sc);

} // namespace mu2e

#endif /* RecoDataProducts_TrackSummary_hh */
