#ifndef CRVDQM_inc_CRVDQMRun1_hh
#define CRVDQM_inc_CRVDQMRun1_hh
// Frozen Run 1 CRV layout for the CRV DQM clients: channel numbering, the
// supported detector configurations with their sectors and counter envelopes,
// and the selection constants that decide what the timing histograms contain.
//
// The CRV DQM histogram set is built from these constants and nothing else, so
// it is the same in every job of a configuration and merges across run
// periods. Every value here is under the clients' kBinningVersion: changing
// one is a version bump.
//
// Configurations (Offline/DQMHelpers/README.md has the sources):
//   run1a      Run 1A/1B, geometry run1a_v01: sectors T1 (17 modules), T2 (3);
//              the full CRVId address space, 18 ROCs / 432 FEB ports.
//   extracted  KPP, geometries extracted_v02..v04: EX, T1, T2 and muon taggers
//              M1-M8 (v04). Envelope is the union over v02-v04. Read out by
//              2 ROCs: ROC 1 fully populated (24 ports), ROC 2 on ports 2-7,
//              so global FEB IDs 0-23 and 25-30 (24 is the empty port 1).
//              ROC 2 is physically ROC 4 on DTC link 3: the unpacker's
//              useROC4asROC2 renumbers its digis, but a status block keeps
//              the physical link, so links 0-3 of DTC 0 are booked and the
//              second ROC's status is on link 3 (links 1 and 2 stay empty).
//
// A client books its FEB-, channel- and link-indexed histograms for one
// configuration (CRVDigiDQM(hists, configuration)); the layout's nROC / nFebs /
// nLinks / nOfflineChannels size those axes through the helpers at the bottom.
// The numbering itself (globalFebId etc.) is the same for every configuration.
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
constexpr int kNFebs = kNROC * kNFebPerROC;               //432
constexpr int kNOnlineChannels = kNFebs * kNChanPerFEB;   //27648
constexpr int kNOfflineChannels = static_cast<int>(CRVId::nChannels);
constexpr int kNChanPerBar = static_cast<int>(CRVId::nChanPerBar);
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
constexpr int globalFebId(int roc, int feb) { return (roc - 1) * kNFebPerROC + feb - 1; }
constexpr int rocChannel(int feb, int febChannel) { return (feb - 1) * kNChanPerFEB + febChannel; }
constexpr int onlineChannel(int roc, int feb, int febChannel)
{
  return globalFebId(roc, feb) * kNChanPerFEB + febChannel;
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
  // Readout booked: ROCs (dense from ROC 1, as the digis number them), FEB
  // ports and DTC links. nFebs is the dense global FEB ID range 0..nFebs-1 and
  // need not be a multiple of kNFebPerROC when the last ROC is partly
  // populated. nLinks is the physical link range 0..nLinks-1
  // (dtcId*nROCPerDTC + linkId, as the status blocks carry it), which can
  // exceed nROC when the ROCs are not on consecutive links.
  int nROC;
  int nFebs;
  int nLinks;
  // Offline channels booked (bar index x kNChanPerBar + SiPM, dense per
  // geometry): those of the configuration's largest geometry.
  int nOfflineChannels;
};

constexpr Layout kLayouts[kNConfigurations] = {
    {"run1a", 2, {"T1", "T2"}, {-6904., 4736., 3962.}, {-904., 4844., 20610.},
     kNROC, kNFebs, kNLinks, 1280 * kNChanPerBar},  //run1a_v01: 1280 bars
    {"extracted", 11,
     {"EX", "T1", "T2", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"},
     {-6904., 4227., 20975.}, {-904., 4719., 24716.},
     2, kNFebPerROC + 7, 4, 576 * kNChanPerBar},  //extracted_v04: 576 bars (v02, v03: 512)
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

// Readout-indexed axes of a configuration. Counts() axes are centred on the
// integer; the *Edges axes run 0..n for the 2D maps that index by bin.
constexpr int nROC(int configuration) { return kLayouts[configuration].nROC; }
constexpr int nFebs(int configuration) { return kLayouts[configuration].nFebs; }
constexpr int nLinks(int configuration) { return kLayouts[configuration].nLinks; }
constexpr int nOnlineChannels(int configuration)
{
  return nFebs(configuration) * kNChanPerFEB;
}
constexpr DQMAxis globalFebAxis(int configuration)
{
  return DQMAxis::Counts(0, nFebs(configuration) - 1);
}
constexpr DQMAxis globalFebEdgesAxis(int configuration)
{
  return DQMAxis(nFebs(configuration), 0., nFebs(configuration));
}
constexpr DQMAxis onlineChannelAxis(int configuration)
{
  return DQMAxis::Counts(0, nOnlineChannels(configuration) - 1);
}
constexpr DQMAxis fpgaPairAxis(int configuration)
{
  return DQMAxis::Counts(0, nFebs(configuration) * kNFpgaPairs - 1);
}
constexpr DQMAxis linkAxis(int configuration)
{
  return DQMAxis::Counts(0, nLinks(configuration) - 1);
}
// Offline-channel axis of a configuration (crvDigisPerChannel, crvPEsVsChannel).
constexpr int nOfflineChannels(int configuration)
{
  return kLayouts[configuration].nOfflineChannels;
}
constexpr DQMAxis offlineChannelAxis(int configuration)
{
  return DQMAxis::Counts(0, nOfflineChannels(configuration) - 1);
}
// A readout address inside the configuration's booked range: filling code
// uses this, not onlineIdInRange, to decide whether a digi has a bin.
constexpr bool febInConfiguration(int configuration, int globalFeb)
{
  return globalFeb >= 0 && globalFeb < nFebs(configuration);
}
constexpr bool linkInConfiguration(int configuration, int link)
{
  return link >= 0 && link < nLinks(configuration);
}
constexpr bool offlineChannelInConfiguration(int configuration, int channel)
{
  return channel >= 0 && channel < nOfflineChannels(configuration);
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
