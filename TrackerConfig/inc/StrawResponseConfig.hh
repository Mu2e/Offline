#ifndef TrackerConditions_StrawResponseConfig_hh
#define TrackerConditions_StrawResponseConfig_hh
//
// Initialize StrawResponse from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct StrawResponseConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence<double>  eDep {
      Name("eDep"), Comment(" KeV ")};
    fhicl::Sequence<double> halfPropVelocity {
      Name("halfPropVelocity"), Comment(" mm/ns ")};
    fhicl::Atom<double> centralWirePos {
      Name("centralWirePos"), Comment(" mm ")};
    fhicl::Sequence<double> tdCentralRes {
      Name("tdCentralRes"), Comment(" tdCentralRes ")};
    fhicl::Sequence<double> tdResSlope {
      Name("tdResSlope"), Comment(" tdResSlope ")};
    fhicl::Sequence<double>  totDriftTime {
      Name("totDriftTime"), Comment(" totDriftTime ")};
    fhicl::Atom<bool>  useDriftErrorCalibration {
      Name("useDriftErrorCalibration"), Comment(" useDriftErrorCalibration ")};
    fhicl::Sequence<double> driftErrorParameters {
      Name("driftErrorParameters"), Comment(" driftErrorParameters ")};
    fhicl::Atom<bool>  useParameterizedDriftErrors {
      Name("useParameterizedDriftErrors"), Comment(" use errors calculated from formula instead of numbers from fcl ")};
    fhicl::Atom<int> parameterizedDriftBins {
      Name("parameterizedDriftBins"), Comment(" number of bins for calculating error and drift offset ")};
    fhicl::Atom<double> parameterizedDriftSigma {
      Name("parameterizedDriftSigma"), Comment(" sigma for calculating drift error and offset ")};
    fhicl::Atom<double> parameterizedDriftTau {
      Name("parameterizedDriftTau"), Comment(" tau for calculating drift error and offset ")};

    fhicl::Atom<double> wireLengthBuffer {
      Name("wireLengthBuffer"), Comment(" wireLengthBuffer ")};
    fhicl::Atom<double> strawLengthFactor {
      Name("strawLengthFactor"), Comment(" strawLengthFactor ")};
    fhicl::Atom<double> errorFactor {
      Name("errorFactor"), Comment(" errorFactor ")};
    fhicl::Atom<bool> useNonLinearDrift {
      Name("useNonLinearDrift"), Comment(" useNonLinearDrift ")};
    fhicl::Atom<double> linearDriftVelocity {
      Name("linearDriftVelocity"), Comment(" mm/ns, only used if nonlindrift= ")};
    fhicl::Atom<double> minDriftRadiusResolution {
      Name("minDriftRadiusResolution"), Comment(" mm ")};
    fhicl::Atom<double> maxDriftRadiusResolution {
      Name("maxDriftRadiusResolution"), Comment(" mm ")};
    fhicl::Atom<double> driftRadiusResolutionRadius {
      Name("driftRadiusResolutionRadius"), Comment(" mm ")};
    fhicl::Atom<double> minT0DOCA {
      Name("minT0DOCA"), Comment("FIXME should be moved to a reconstruction configuration ")};
    fhicl::Atom<double> t0shift {
      Name("t0shift"), Comment("FIXME should be average slewing?")};
    fhicl::Sequence<double> peakMinusPedestalEnergyScale {
      Name("peakMinusPedestalEnergyScale"), Comment(" fudge factor for peak minus pedestal energy method ")};
    fhicl::Sequence<double> timeOffsetPanel {
      Name("timeOffsetPanel"), Comment(" timeOffsetPanel ")};
    fhicl::Sequence<double> timeOffsetStrawHV {
      Name("timeOffsetStrawHV"), Comment(" timeOffsetStrawHV ")};
    fhicl::Sequence<double> timeOffsetStrawCal {
      Name("timeOffsetStrawCal"), Comment(" timeOffsetStrawCal ")};

    fhicl::OptionalAtom<double> electronicsTimeDelay {
      Name("electronicsTimeDelay"), Comment("electronicsTimeDelay")};
    fhicl::OptionalAtom<double> gasGain {
      Name("gasGain"), Comment("gasGain")};
    fhicl::OptionalAtom<double> thresholdAnalogNoise {
      Name("thresholdAnalogNoise"), Comment("thresholdAnalogNoise")};
    fhicl::OptionalAtom<double> adcAnalogNoise {
      Name("adcAnalogNoise"), Comment("adcAnalogNoise")};
    fhicl::OptionalAtom<double> defaultThresholddVdI {
      Name("defaultThresholddVdI"), Comment("defaultThresholddVdI")};
    fhicl::OptionalAtom<double> defaultAdcdVdI {
      Name("defaultAdcdVdI"), Comment("defaultAdcdVdI")};
    fhicl::OptionalAtom<double> saturationVoltage {
      Name("saturationVoltage"), Comment("saturationVoltage")};
    fhicl::OptionalAtom<double> ADCPedestal {
      Name("ADCPedestal"), Comment("ADCPedestal")};

  };

}

#endif
