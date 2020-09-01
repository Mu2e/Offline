#ifndef TrackerConditions_StrawDriftConfig_hh
#define TrackerConditions_StrawDriftConfig_hh
//
// Initialize StrawDrift from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct StrawDriftConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Atom<double> wireVoltage {
      Name("wireVoltage"), Comment("wire voltage in V")};
    fhicl::Atom<int> phiBins{
      Name("phiBins"), Comment("number of bins in phi for drift model")};
    fhicl::Atom<double> deltaDistance{
      Name("deltaDistance"), Comment("Size of bins in distance for drift D2T model")};
    fhicl::Atom<double> deltaTime{
      Name("deltaTime"), Comment("Size of bins in time for drift T2D model")};
    fhicl::Atom<int> driftIntegrationBins{
      Name("driftIntegrationBins"), 
	Comment("number of integrations steps for drift model")};
    fhicl::Sequence<double> kVcm{
      Name("kVcm"), Comment("drift model field in KV/cm") };
    fhicl::Sequence<double> cmus{
      Name("cmus"), Comment("drift model speed in cm/us") };
  };

}

#endif
