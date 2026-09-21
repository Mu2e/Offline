#ifndef CRVDQM_inc_CRVDQMRun1_hh
#define CRVDQM_inc_CRVDQMRun1_hh
// Frozen Run 1 CRV layout for the CRV DQM clients: channel numbering, the
// supported detector configurations with their sectors and counter envelopes,
// and the selection constants that decide what the timing histograms contain.
//
// The CRV DQM histogram set is built from these constants and nothing else, so
// it is the same in every job and merges across run periods. Every value here
// is under the clients' kBinningVersion: changing one is a version bump.
//
// Configurations (Offline/DQMHelpers/README.md has the sources):
//   run1a      Run 1A/1B, geometry run1a_v01: sectors T1 (17 modules), T2 (3)
//   extracted  KPP, geometries extracted_v02..v04: EX, T1, T2 and muon taggers
//              M1-M8 (v04). Envelope is the union over v02-v04.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMAxis.hh"
#include "Offline/DataProducts/inc/CRVId.hh"

#include <string>

namespace mu2e {
namespace CRVDQMRun1 {

// Dense CRVId numbering. ROC and FEB are 1-based on the wire, FEB channel is 0-based.
constexpr int kNROC = static_cast<int>(CRVId::nROC);
constexpr int kNFebPerROC = static_cast<int>(CRVId::nFEBPerROC);
constexpr int kNChanPerFEB = static_cast<int>(CRVId::nChanPerFEB);
constexpr int kNChanPerROC = kNFebPerROC * kNChanPerFEB;
constexpr int kNFebPorts = kNROC * kNFebPerROC;               //432
constexpr int kNOnlineChannels = kNFebPorts * kNChanPerFEB;   //27648
constexpr int kNOfflineChannels = static_cast<int>(CRVId::nChannels);
constexpr int kNFPGAPerFEB = static_cast<int>(CRVId::nFPGAPerFEB);
constexpr int kNChanPerFPGA = static_cast<int>(CRVId::nChanPerFPGA);
// Unordered FPGA pairs on one FEB, same FPGA included: 10.
constexpr int kNFpgaPairs = kNFPGAPerFEB * (kNFPGAPerFEB + 1) / 2;
// DTC links: dtcId * nROCPerDTC + linkId.
constexpr int kNLinksPerDTC = static_cast<int>(CRVId::nROCPerDTC);
constexpr int kNLinks = kNROC;

constexpr bool onlineIdInRange(int roc, int feb, int febChannel)
{
  return roc >= 1 && roc <= kNROC && feb >= 1 && feb <= kNFebPerROC &&
         febChannel >= 0 && febChannel < kNChanPerFEB;
}
constexpr int febPort(int roc, int feb) { return (roc - 1) * kNFebPerROC + feb - 1; }
constexpr int rocChannel(int feb, int febChannel) { return (feb - 1) * kNChanPerFEB + febChannel; }
constexpr int onlineChannel(int roc, int feb, int febChannel)
{
  return febPort(roc, feb) * kNChanPerFEB + febChannel;
}
// fpgaA <= fpgaB
constexpr int fpgaPairIndex(int fpgaA, int fpgaB)
{
  return fpgaA * kNFPGAPerFEB - fpgaA * (fpgaA - 1) / 2 + (fpgaB - fpgaA);
}
constexpr bool linkInRange(int dtcId, int linkId)
{
  return linkId >= 0 && linkId < kNLinksPerDTC && dtcId >= 0 &&
         dtcId * kNLinksPerDTC + linkId < kNLinks;
}
constexpr int globalLink(int dtcId, int linkId) { return dtcId * kNLinksPerDTC + linkId; }

// Constant-fraction timing of a digi waveform.
constexpr double kCFFraction = 0.20;
constexpr int kCFMinAmplitude = 10;       //ADC above the first sample

// Partner-FEB timing: a hit counts at this amplitude or above (~11 PE at the
// KPP calibration), and a local group is at least kDtMinLayers of a module
// group's 4 layers with consecutive hits no more than kDtCoincWindow apart.
constexpr int kDtMinAmplitude = 200;      //ADC
constexpr double kDtCoincWindow = 20.;    //ns
constexpr int kDtMinLayers = 3;

enum Configuration { kRun1a = 0, kExtracted, kNConfigurations };
constexpr int kMaxSectors = 11;

struct Layout {
  const char* name;
  int nSectors;
  const char* sectors[kMaxSectors];
  // Counter envelope in Mu2e coordinates [mm], measured from each geometry.
  double lo[3];
  double hi[3];
};

constexpr Layout kLayouts[kNConfigurations] = {
    {"run1a", 2, {"T1", "T2"}, {-6904., 4736., 3962.}, {-904., 4844., 20610.}},
    {"extracted", 11,
     {"EX", "T1", "T2", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"},
     {-6904., 4227., 20975.}, {-904., 4719., 24716.}},
};

// Geometry crs.name -> configuration. extracted_v01..v03 all say "extracted".
struct GeometryAlias {
  const char* crsName;
  Configuration configuration;
};
constexpr GeometryAlias kGeometries[] = {
    {"run1a_v01", kRun1a},
    {"extracted", kExtracted},
    {"extracted_v04", kExtracted},
};

// Cluster position axes: the envelope plus 5% of its span on each side.
constexpr int kNPosBins = 100;
constexpr double kPosMargin = 0.05;
constexpr DQMAxis positionAxis(int configuration, int coordinate)
{
  const Layout& l = kLayouts[configuration];
  const double pad = kPosMargin * (l.hi[coordinate] - l.lo[coordinate]);
  return DQMAxis(kNPosBins, l.lo[coordinate] - pad, l.hi[coordinate] + pad);
}

const char* configurationName(int configuration);
// -1 for a geometry this layout does not know: the caller must throw.
int configurationFor(const std::string& crsName);
// Index of `sector` ("EX" or the geometry's "CRV_EX") in the configuration's list, or -1.
int sectorIndex(int configuration, const std::string& sector);
// Histogram name suffix, e.g. "_extracted_EX".
std::string sectorTag(int configuration, int sector);

} // namespace CRVDQMRun1
} // namespace mu2e

#endif /* CRVDQM_inc_CRVDQMRun1_hh */
