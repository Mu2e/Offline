#ifndef TrackerConditions_StrawPhysicsConfig_hh
#define TrackerConditions_StrawPhysicsConfig_hh
//
// Initialize StrawPhysics from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct StrawPhysicsConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Atom<double>  meanFreePath {
      Name("meanFreePath"), Comment(" mm, average distance between ionizations for a MIP in STP Ar (Blum etal, table 1.1)")};
    fhicl::Atom<double>  ionizedElectronKE {
      Name("ionizedElectronKE"), Comment("kinetic energy of electron (MeV)")};
    fhicl::Atom<double> electronCharge {
      Name("electronCharge"), Comment("e, pC")};
    fhicl::Atom<double> gasGain {
      Name("gasGain"), Comment("gas gain")};
    fhicl::Atom<double> polyaA {
      Name("polyaA"), Comment("A = 1/f = theta + 1.  A=1 -> exponential, A=infinity->delta-function")};
    fhicl::Atom<double> gainRMSSlope {
      Name("gainRMSSlope"), Comment("slope of relative gain sigma on 1/sqrt(n)")};
    fhicl::Atom<double> nGainGauss {
      Name("nGainGauss"), Comment("number of electrons/cluster to switch to a Gaussian model of the gain fluctuations")};
    fhicl::Atom<double> propagationVelocity {
      Name("propagationVelocity"), Comment("propagation velocity")};
    fhicl::Sequence<double> clusterDriftPolynomial {
      Name("clusterDriftPolynomial"), Comment("linear term has units nanoseconds/mm")};
    fhicl::Sequence<double> driftTimeVariance {
      Name("driftTimeVariance"), Comment("Drift time variance linear dependence on drift time (ns)")};
    fhicl::Atom<bool> useNonLinearDrift {
      Name("useNonLinearDrift"), Comment("switch to turn on/off non linear drift for diagnosis")};
    fhicl::Atom<double> bFieldOverride {
      Name("bFieldOverride"), Comment("")};

    fhicl::Sequence<double> probPerCharge {
      Name("probPerCharge"), Comment("Blum, table 1.4")};
    fhicl::Sequence<double> ionizationEnergyTable {
      Name("ionizationEnergyTable"), Comment("")};



  };

}

#endif
