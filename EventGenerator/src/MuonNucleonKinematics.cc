// Kinematics helpers for bound mu- + nucleon annihilation, see the header for details
//
// Michael MacKenzie, 2026

#include "Offline/EventGenerator/inc/MuonNucleonKinematics.hh"

#include <algorithm>
#include <cmath>

#include "cetlib_except/exception.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  //================================================================
  MuonNucleonKinematics::Model MuonNucleonKinematics::parseModel(const std::string& name) {
    if(name == "free")     return Model::free;
    if(name == "fermiGas") return Model::fermiGas;
    if(name == "nucleus")  return Model::nucleus;
    throw cet::exception("BADCONFIG") << "MuonNucleonKinematics: nuclear model \"" << name
                                      << "\" not defined, only \"free\", \"fermiGas\", or \"nucleus\" are allowed\n";
  }

  //================================================================
  MuonNucleonKinematics::NuclearParameters MuonNucleonKinematics::parseConfig(const NuclearConfig& conf, const std::string& caller) {
    NuclearParameters par{parseModel(conf.nuclearModel()), conf.nuclearMaterial()};
    const bool hasKF = conf.fermiMomentum(par.fermiMomentum);
    const bool hasRemoval = conf.removalEnergy(par.removalEnergy);
    const bool needKF = par.model == Model::fermiGas;
    const bool needRemoval = par.model != Model::free;
    if(hasKF != needKF) {
      throw cet::exception("BADCONFIG") << caller << ": fermiMomentum must " << ((needKF) ? "" : "not ")
                                        << "be set for nuclearModel \"" << conf.nuclearModel() << "\"\n";
    }
    if(hasRemoval != needRemoval) {
      throw cet::exception("BADCONFIG") << caller << ": removalEnergy must " << ((needRemoval) ? "" : "not ")
                                        << "be set for nuclearModel \"" << conf.nuclearModel() << "\"\n";
    }
    if(par.fermiMomentum < 0.) throw cet::exception("BADCONFIG") << caller << ": negative Fermi momentum\n";
    if(par.removalEnergy < 0.) throw cet::exception("BADCONFIG") << caller << ": negative removal energy\n";
    return par;
  }

  //================================================================
  void MuonNucleonKinematics::checkMaterial(const NuclearParameters& par, const std::string& material, const std::string& caller) {
    if(par.material != material) {
      throw cet::exception("BADCONFIG") << caller << ": nuclear parameters configured for \"" << par.material
                                        << "\" but the stopping target material is \"" << material << "\"\n";
    }
  }

  //================================================================
  MuonNucleonKinematics::MuonNucleonKinematics(CLHEP::HepRandomEngine& eng,
                                               double muonEnergy, double nucleonMass, double nucleusMass,
                                               const NuclearParameters& par) :
    _muonEnergy(muonEnergy),
    _nucleonMass(nucleonMass),
    _nucleusMass(nucleusMass),
    _residualMass(nucleusMass - nucleonMass + par.removalEnergy),
    _model(par.model),
    _fermiMomentum(par.fermiMomentum),
    _removalEnergy(par.removalEnergy),
    _randFlat(eng),
    _randomUnitSphere(eng)
  {}

  //================================================================
  double MuonNucleonKinematics::maxInvariantMass() const {
    return _muonEnergy + _nucleonMass - ((_model == Model::free) ? 0. : _removalEnergy);
  }

  //================================================================
  void MuonNucleonKinematics::generate(const std::vector<double>& masses, std::vector<CLHEP::HepLorentzVector>& momenta,
                                       double* nucleonMomentum) {
    double massSum = 0.;
    for(const double m : masses) massSum += m;
    if(massSum >= maxInvariantMass()) {
      throw cet::exception("BADCONFIG") << "MuonNucleonKinematics: final state mass sum " << massSum
                                        << " MeV exceeds the maximum available invariant mass " << maxInvariantMass() << " MeV\n";
    }

    switch(_model) {
    case Model::free: {
      if(nucleonMomentum) *nucleonMomentum = 0.;
      phaseSpace(CLHEP::HepLorentzVector(0., 0., 0., _muonEnergy + _nucleonMass), masses, momenta);
      return;
    }
    case Model::nucleus: {
      // Residual nucleus as an extra final state body
      std::vector<double> all(masses);
      all.push_back(_residualMass);
      phaseSpace(CLHEP::HepLorentzVector(0., 0., 0., _muonEnergy + _nucleusMass), all, momenta);
      if(nucleonMomentum) *nucleonMomentum = momenta.back().vect().mag();
      momenta.pop_back();
      return;
    }
    case Model::fermiGas: {
      // Uniform in the Fermi sphere: p^2 dp --> p = k_F * u^(1/3). Reject states below threshold.
      constexpr int maxTries = 1000000;
      for(int tries = 0; tries < maxTries; ++tries) {
        const double p = _fermiMomentum * std::cbrt(_randFlat.fire());
        const CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(p);
        const double recoil = p*p/(2.*_residualMass);
        const CLHEP::HepLorentzVector total(p3, _muonEnergy + _nucleonMass - _removalEnergy - recoil);
        if(total.m() <= massSum) continue;
        if(nucleonMomentum) *nucleonMomentum = p;
        phaseSpace(total, masses, momenta);
        return;
      }
      throw cet::exception("RUNTIME") << "MuonNucleonKinematics: failed to sample an initial state above threshold\n";
    }
    }
  }

  //================================================================
  double MuonNucleonKinematics::pdk(double a, double b, double c) {
    const double x = (a - b - c)*(a + b + c)*(a - b + c)*(a + b - c);
    return (x > 0.) ? std::sqrt(x)/(2.*a) : 0.;
  }

  //================================================================
  // Maximum phase space weight for a given kinetic energy release, following TGenPhaseSpace
  double MuonNucleonKinematics::maxWeight(double teCM, const std::vector<double>& masses) const {
    double emmax = teCM + masses[0];
    double emmin = 0.;
    double wtMax = 1.;
    for(size_t n = 1; n < masses.size(); ++n) {
      emmin += masses[n-1];
      emmax += masses[n];
      wtMax *= pdk(emmax, emmin, masses[n]);
    }
    return wtMax;
  }

  //================================================================
  void MuonNucleonKinematics::phaseSpace(const CLHEP::HepLorentzVector& total, const std::vector<double>& masses,
                                         std::vector<CLHEP::HepLorentzVector>& momenta) {
    const size_t n = masses.size();
    if(n < 2) throw cet::exception("BADCONFIG") << "MuonNucleonKinematics: phase space needs at least two particles\n";
    double massSum = 0.;
    for(const double m : masses) massSum += m;
    const double teCM = total.m() - massSum;
    if(teCM <= 0.) throw cet::exception("RUNTIME") << "MuonNucleonKinematics: invariant mass " << total.m()
                                                   << " below the final state mass sum " << massSum << "\n";

    std::vector<double> invMas(n), pd(n), rno(n);
    momenta.assign(n, CLHEP::HepLorentzVector());

    // Accept-reject on the phase space weight at fixed invariant mass, so the initial state
    // distribution follows the nuclear model unbiased by the available phase space volume
    const double wtMax = maxWeight(teCM, masses);
    constexpr int maxTries = 1000000;
    int tries = 0;
    while(true) {
      if(++tries > maxTries) throw cet::exception("RUNTIME") << "MuonNucleonKinematics: failed to sample phase space\n";
      rno[0] = 0.;
      rno[n-1] = 1.;
      for(size_t i = 1; i + 1 < n; ++i) rno[i] = _randFlat.fire();
      std::sort(rno.begin() + 1, rno.end() - 1);
      double sum = 0.;
      for(size_t i = 0; i < n; ++i) {
        sum += masses[i];
        invMas[i] = rno[i]*teCM + sum;
      }
      double wt = 1.;
      for(size_t i = 0; i + 1 < n; ++i) {
        pd[i] = pdk(invMas[i+1], invMas[i], masses[i+1]);
        wt *= pd[i];
      }
      if(_randFlat.fire()*wtMax <= wt) break;
    }

    // Build the momenta in successive rest frames
    momenta[0].setVectM(CLHEP::Hep3Vector(0.,  pd[0], 0.), masses[0]);
    momenta[1].setVectM(CLHEP::Hep3Vector(0., -pd[0], 0.), masses[1]);
    size_t i = 1;
    while(true) {
      // Random isotropic orientation of the subsystem
      const double angZ = std::acos(2.*_randFlat.fire() - 1.);
      const double angY = CLHEP::twopi*_randFlat.fire();
      for(size_t j = 0; j <= i; ++j) {
        momenta[j].rotateZ(angZ);
        momenta[j].rotateY(angY);
      }
      if(i == n - 1) break;
      const double beta = pd[i] / std::sqrt(pd[i]*pd[i] + invMas[i]*invMas[i]);
      for(size_t j = 0; j <= i; ++j) momenta[j].boost(CLHEP::Hep3Vector(0., beta, 0.));
      ++i;
      momenta[i].setVectM(CLHEP::Hep3Vector(0., -pd[i-1], 0.), masses[i]);
    }

    // Boost to the lab frame
    const CLHEP::Hep3Vector boost = total.boostVector();
    for(auto& mom : momenta) mom.boost(boost);
  }

  //================================================================
  void MuonNucleonKinematics::twoBodyDecay(const CLHEP::HepLorentzVector& parent, double m1, double m2,
                                           CLHEP::HepLorentzVector& p1, CLHEP::HepLorentzVector& p2) {
    const double p = pdk(parent.m(), m1, m2);
    const CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(p);
    p1.setVectM( p3, m1);
    p2.setVectM(-p3, m2);
    const CLHEP::Hep3Vector boost = parent.boostVector();
    p1.boost(boost);
    p2.boost(boost);
  }
}
