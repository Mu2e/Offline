#ifndef EventGenerator_MuonNucleonKinematics_hh
#define EventGenerator_MuonNucleonKinematics_hh
//
// Kinematics helpers for the annihilation of a bound mu- with a bound nucleon in a muonic atom,
// used by the baryon number violating muon capture generator tools.
//
// The muon is at rest with its bound 1s total energy. Three models of the nucleon are available:
//  - fermiGas: Fermi gas spectator model. The nucleon momentum is uniform within the Fermi sphere,
//    its energy is m_N - removalEnergy, and the residual (A-1) nucleus recoils with the opposite
//    momentum and p^2/(2 M_{A-1}) of kinetic energy. Initial states below the final state mass
//    threshold are rejected.
//  - nucleus: the residual nucleus (mass M_A - m_N + removalEnergy, i.e. removalEnergy is the
//    nucleon separation energy plus any excitation) is included as an extra body in the final
//    state phase space with total energy E_mu + M_A, as in arXiv:2407.03450.
//  - free: the nucleon is at rest with no binding (muonic hydrogen-like).
//
// Final state: n-body Lorentz-invariant phase space at fixed invariant mass, following the
// Raubold-Lynch method (F. James, CERN 68-15) as implemented in ROOT's TGenPhaseSpace, but sampled
// with the art random engine.
//
// Michael MacKenzie, 2026

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

namespace mu2e {

  class MuonNucleonKinematics {
  public:
    enum class Model { free, fermiGas, nucleus };
    static Model parseModel(const std::string& name);

    // nucleusMass is the mass of the initial nucleus M_A, used for the residual nucleus
    MuonNucleonKinematics(CLHEP::HepRandomEngine& eng,
                          double muonEnergy, double nucleonMass, double nucleusMass,
                          Model model, double fermiMomentum, double removalEnergy);

    // Sample the final state of mu- + nucleon --> particles with the given masses, in the lab frame.
    // If nucleonMomentum is given, it is filled with the momentum magnitude of the struck nucleon
    // (the residual nucleus recoil momentum).
    void generate(const std::vector<double>& masses, std::vector<CLHEP::HepLorentzVector>& momenta,
                  double* nucleonMomentum = nullptr);

    // Maximum invariant mass available to the annihilation products
    double maxInvariantMass() const;

    // Unweighted n-body phase space (n >= 2) at the invariant mass of total, boosted to its frame
    void phaseSpace(const CLHEP::HepLorentzVector& total, const std::vector<double>& masses,
                    std::vector<CLHEP::HepLorentzVector>& momenta);

    // Isotropic two-body decay of parent --> 1 + 2 in the parent rest frame, boosted to the lab frame
    void twoBodyDecay(const CLHEP::HepLorentzVector& parent, double m1, double m2,
                      CLHEP::HepLorentzVector& p1, CLHEP::HepLorentzVector& p2);

    // Two-body decay momentum of a --> b + c
    static double pdk(double a, double b, double c);

  private:
    double maxWeight(double teCM, const std::vector<double>& masses) const;

    double _muonEnergy;
    double _nucleonMass;
    double _nucleusMass;
    double _residualMass;
    Model  _model;
    double _fermiMomentum;
    double _removalEnergy;

    CLHEP::RandFlat  _randFlat;
    RandomUnitSphere _randomUnitSphere;
  };
}

#endif
