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
    fhicl::Atom<double> eBins {
      Name("eBins"), Comment("Number of energy bins")};
    fhicl::Atom<double> eBinWidth {
      Name("eBinWidth"), Comment("Width of energy bins (KeV)")};
    fhicl::Sequence<double> ehalfPVScale {
      Name("eHalfPVScale"), Comment(" fraction of nominal ")};
    fhicl::Atom<double> centralWirePos {
      Name("centralWirePos"), Comment(" mm ")};
    fhicl::Sequence<double> tdCentralRes {
      Name("tdCentralRes"), Comment(" tdCentralRes ")};
    fhicl::Sequence<double> tdResSlope {
      Name("tdResSlope"), Comment(" tdResSlope ")};
    fhicl::Sequence<double> strawHalfPropVelocity {
      Name("strawHalfPropVelocity"), Comment(" mm/ns ")};
    fhicl::Atom<double> defaultHalfPropVelocity {
      Name("defaultHalfPropVelocity"), Comment(" mm/ns ")};
    fhicl::Atom<bool> truncateLongitudinal {
      Name("truncateLongitudinal"), Comment("truncate reco longitudinal at straw end")};
    fhicl::Atom<bool> rmsLongErrors {
      Name("rmsLongErrors"), Comment("Use errors tuned from residual profile rms")};
    fhicl::Atom<int> totTBins {
      Name("totTBins"), Comment("TOT drift time t bins")};
    fhicl::Atom<double> totTBinWidth {
      Name("totTBinWidth"), Comment("TOT drift time t bin width")};
    fhicl::Atom<int> totEBins {
      Name("totEBins"), Comment("TOT drift time e bins")};
    fhicl::Atom<double> totEBinWidth {
      Name("totEBinWidth"), Comment("TOT drift time e bin width")};
    fhicl::Sequence<double>  totDriftTime {
      Name("totDriftTime"), Comment(" totDriftTime ")};
    fhicl::Sequence<double>  totDriftError {
      Name("totDriftError"), Comment(" totDriftError ")};
    fhicl::Sequence<double> llDriftTimeOffBins {
      Name("llDriftTimeOffBins"), Comment("Drift time Offset Bin edges for likelihood fit (mm)")};
    fhicl::Sequence<double> llDriftTimeOffset {
      Name("llDriftTimeOffset"), Comment("Drift time offset for likelihood fit (ns)")};
    fhicl::Sequence<double> llDriftTimeRMSBins {
      Name("llDriftTimeRMSBins"), Comment("Drift time RMS Bin edges for likelihood fit (mm)")};
    fhicl::Sequence<double> llDriftTimeRMS {
      Name("llDriftTimeRMS"), Comment("Drift time RMS for likelihood fit (ns)")};
    fhicl::Sequence<double> driftOffBins {
      Name("driftOffBins"), Comment("Drift Offset Bin edges (mm)")};
    fhicl::Sequence<double> driftOffset {
      Name("driftOffset"), Comment("Drift Offset (mm)")};
    fhicl::Sequence<double> driftRMSBins {
      Name("driftRMSBins"), Comment("Drift RMS Bin edges (mm)")};
    fhicl::Sequence<double> signedDriftRMS {
      Name("signedDriftRMS"), Comment("Signed Drift RMS (mm)")};
    fhicl::Sequence<double> unsignedDriftRMS {
      Name("unsignedDriftRMS"), Comment("Unsigned Drift RMS (mm)")};
    fhicl::Atom<double> dRdTScale {
      Name("dRdTScale"), Comment("Scale factor for dRdT")};

    fhicl::Atom<bool> driftIgnorePhi {
      Name("driftIgnorePhi"), Comment("Ignore phi for no field reco")};

    fhicl::Atom<double> wireLengthBuffer {
      Name("wireLengthBuffer"), Comment(" wireLengthBuffer ")};
    fhicl::Atom<double> strawLengthFactor {
      Name("strawLengthFactor"), Comment(" strawLengthFactor ")};
    fhicl::Atom<double> errorFactor {
      Name("errorFactor"), Comment("Error scaling for longitudinal reco outside straw length")};
    fhicl::Atom<bool> useNonLinearDrift {
      Name("useNonLinearDrift"), Comment(" useNonLinearDrift ")};
    fhicl::Atom<double> linearDriftVelocity {
      Name("linearDriftVelocity"), Comment(" mm/ns, only used if nonlindrift= ")};
    fhicl::Atom<double> defaultPeakMinusPedestalEnergyScale {
      Name("defaultPeakMinusPedestalEnergyScale"), Comment("default constant value for pmp energy method calibration")};
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
