#ifndef TrackerConditions_StrawElectronicsConfig_hh
#define TrackerConditions_StrawElectronicsConfig_hh
//
// Initialize StrawElectronics from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct StrawElectronicsConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 

    fhicl::Atom<double> defaultThresholddVdI {
      Name("defaultThresholddVdI"), Comment("defaultThresholddVdI")};
    fhicl::Sequence<double>  thresholddVdI {
      Name("thresholddVdI"), Comment("threshold dVdI")};
    fhicl::Atom<double> defaultAdcdVdI {
      Name("defaultAdcdVdI"), Comment("default Adc dVdI")};
    fhicl::Sequence<double>  adcdVdI {
      Name("adcdVdI"), Comment("Adc dVdI, mVolt/uAmps (transimpedance gain)")};

    fhicl::Atom<double> deadTimeAnalog {
      Name("deadTimeAnalog"), Comment("nsec dead after threshold crossing (pulse baseline restoration time)")};
    fhicl::Atom<double> deadTimeDigital {
      Name("deadTimeDigital"), Comment("sec dead after threshold crossing (electronics processing time)")};
    fhicl::Atom<double> saturationVoltage {
      Name("saturationVoltage"), Comment("mVolt")};
    fhicl::Atom<double> defaultDiscriminatorThreshold {
      Name("defaultDiscriminatorThreshold"), Comment("mVolt, post amplification")};
    fhicl::Sequence<double> discriminatorThreshold {
      Name("discriminatorThreshold"), Comment("mVolt, post amplification")};
    fhicl::Atom<double> strawNoise {
      Name("strawNoise"), Comment("mVolt")};
    fhicl::Atom<double> thresholdAnalogNoise {
      Name("thresholdAnalogNoise"), Comment("mVolt")};
    fhicl::Atom<double> adcAnalogNoise {
      Name("adcAnalogNoise"), Comment("")};
    fhicl::Atom<double> ADCLSB {
      Name("ADCLSB"), Comment("mVolt")};
    fhicl::Atom<int> maxADC {
      Name("maxADC"), Comment("maxADC")};
    fhicl::Atom<unsigned> nADCPresamples {
      Name("nADCPresamples"), Comment("nADCPresamples")};
    fhicl::Atom<double> ADCPeriod {
      Name("ADCPeriod"), Comment("nsec")};
    fhicl::Atom<double> ADCOffset {
      Name("ADCOffset"), Comment("nsec")};
    fhicl::Atom<unsigned> maxThreshTimeSeparation {
      Name("maxThreshTimeSeparation"), Comment("ADC clock ticks")};
    fhicl::Atom<unsigned> tCoince {
      Name("tCoince"), Comment("maxing threshold xing pair time separation to create a digi, in number of ADC clock cycles")};
    fhicl::Atom<double> TDCLSB {
      Name("TDCLSB"), Comment("nsec, least-significant bit of TDC")};
    fhicl::Atom<unsigned> maxTDC {
      Name("maxTDC"), Comment("16 bits, maximum TDC value")};
    fhicl::Atom<double>  TOTLSB {
      Name("TOTLSB"), Comment("ns, least-significant bit of TOT")};
    fhicl::Atom<unsigned> maxTOT {
      Name("maxTOT"), Comment("maximum TOT value")};
    fhicl::Atom<double> TDCResolution {
      Name("TDCResolution"), Comment("nsec, tdc resolution (electronics effects only)")};
    fhicl::Atom<double> electronicsTimeDelay {
      Name("electronicsTimeDelay"), Comment("nsec, Absolute time delay in electronics due to firmware signal propagation etc")};
    fhicl::Atom<double> eventWindowMarkerROCJitter {
      Name("eventWindowMarkerROCJitter"), Comment("ps (jitter per panel per microbuncH)")};
    fhicl::Atom<double> flashStart {
      Name("flashStart"), Comment("nsec, flash blanking period")};
    fhicl::Atom<double> flashEnd {
      Name("flashEnd"), Comment("nsec, flash blanking period")};
    fhicl::Atom<double> flashClockSpeed {
      Name("flashClockSpeed"), Comment("nsec")};
    fhicl::Atom<int> responseBins {
      Name("responseBins"), Comment("")};
    fhicl::Atom<double> sampleRate {
      Name("sampleRate"), Comment("ghz")};
    fhicl::Atom<int> saturationSampleFactor {
      Name("saturationSampleFactor"), Comment("")};
    fhicl::Sequence<double> preampPoles {
      Name("preampPoles"), Comment("preampPoles")};
    fhicl::Sequence<double> preampZeros {
      Name("preampZeros"), Comment("preampZeros")};
    fhicl::Sequence<double> ADCPoles {
      Name("ADCPoles"), Comment("ADCPoles")};
    fhicl::Sequence<double> ADCZeros {
      Name("ADCZeros"), Comment("ADCZeros")};
    fhicl::Sequence<double> preampToAdc1Poles {
      Name("preampToAdc1Poles"), Comment("preampToAdc1Poles")};
    fhicl::Sequence<double> preampToAdc1Zeros {
      Name("preampToAdc1Zeros"), Comment("preampToAdc1Zeros")};
    fhicl::Sequence<double> preampToAdc2Poles {
      Name("preampToAdc2Poles"), Comment("preampToAdc2Poles")};
    fhicl::Sequence<double> preampToAdc2Zeros {
      Name("preampToAdc2Zeros"), Comment("preampToAdc2Zeros")};


    fhicl::Sequence<double> wireDistances {
      Name("wireDistances"), Comment("wireDistances")};
    fhicl::Sequence<double> currentMeans {
      Name("currentMeans"), Comment("currentMeans")};
    fhicl::Sequence<double> currentNormalizations {
      Name("currentNormalizations"), Comment("currentNormalizations")};
    fhicl::Sequence<double> currentSigmas {
      Name("currentSigmas"), Comment("currentSigmas")};
    fhicl::Sequence<double> currentT0s {
      Name("currentT0s"), Comment("currentT0s")};
    fhicl::Atom<double> clusterLookbackTime {
      Name("clusterLookbackTime"), Comment("clusterLookbackTime")};
    fhicl::Sequence<double> timeOffsetPanel {
      Name("timeOffsetPanel"), Comment("timeOffsetPanel")};
    fhicl::Sequence<double> timeOffsetStrawHV {
      Name("timeOffsetStrawHV"), Comment("timeOffsetStrawHV")};
    fhicl::Sequence<double> timeOffsetStrawCal {
      Name("timeOffsetStrawCal"), Comment("timeOffsetStrawCal")};
  };

}

#endif
