// Andrei Gaponenko, 2014

#include "Offline/RecoDataProducts/inc/TrackSummary.hh"

#include <cmath>

#include "CLHEP/Matrix/Vector.h"

#include "KinKal/General/Chisq.hh"

namespace mu2e {
  /* TODO: Delete if OK. Also check if can use BTrkLegacy HelixParams to substitute
  TrackSummary::HelixParams::HelixParams(const TrkSimpTraj& ltraj)
    : d0_(ltraj.parameters()->parameter()[d0Index])
    , phi0_(ltraj.parameters()->parameter()[phi0Index])
    , omega_(ltraj.parameters()->parameter()[omegaIndex])
    , z0_(ltraj.parameters()->parameter()[z0Index])
    , tanDip_(ltraj.parameters()->parameter()[tanDipIndex])
    , covariance_(ltraj.parameters()->covariance())
  {}
  */
  double TrackSummary::HelixParams::radius() const {
    // can check for division by zero here
    return 1./omega_;
  }

  double TrackSummary::HelixParams::dOut() const {
    return d0_ + 2.*radius();
  }

  double TrackSummary::HelixParams::wavelength() const {
    static const double pi = 4.*atan(1.); // M_PI is not in the standard
    return 2 * pi * radius() * tanDip_;
  }

  double TrackSummary::fitcon() const {
    KinKal::Chisq cs(chi2(), ndof());
    return cs.probability();
  }

  void TrackSummary::addState(const TrackSummary::TrackStateAtPoint& st) {
    states_.emplace_back(st);
  }

  double TrackSummary::TrackStateAtPoint::momentumError() const {
    using namespace CLHEP;
    Hep3Vector u(momentum_.unit());
    HepVector dir(3);
    for(unsigned i=0; i<3; ++i) {
      dir[i] = u[i];
    }
    return sqrt(momentumCovariance_.similarity(dir));
  }

  double TrackSummary::TrackStateAtPoint::costh() const {
    const double p = momentum_.mag();
    return (p > 0.) ? momentum_.z()/p : 0.;
  }

  std::ostream& operator<<(std::ostream& os, const TrackSummary::TrackStateAtPoint& st) {
    os<<"(flightLength="<<st.flightLength()
      <<", arrivalTime="<<st.arrivalTime()
      <<", pos="<<st.position()
      <<", mom="<<st.momentum()
      <<", momErr="<<st.momentumError()
      <<" )";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const TrackSummary& sum) {
    os<<"TrackSummary(q="<<sum.charge()
      <<", nactive="<<sum.nactive()
      <<", ndof="<<sum.ndof()
      <<", chi2="<<sum.chi2()
      <<", t0="<<sum.t0()<<" +- "<<sum.t0Err()
      <<", flt0="<<sum.flt0()
      <<", states = [";

    for(const auto& i: sum.states()) {
      os<<i<<", ";
    }

    os<<" ] )";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const TrackSummaryCollection& sc) {
    os<<"TrackSummaryCollection: "<<sc.size()<<" entries\n";
    for(const auto& i: sc) {
      os<<"\t"<<i<<"\n";
    }
    return os;
  }

}
