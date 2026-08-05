#pragma once

// File-name convention for the VD resampler pipeline (single source of truth).
//
// Names follow the Mu2e dataset convention <tier>.mu2e.<description>.<config>.<sequencer>.<ext>,
// where the sequencer is <run6>_<src8>: a 6-digit run number and the 8-digit enumeration of the
// data source. Every artifact embeds versionTag, run number and dataSourceTag so several campaigns
// can coexist in one directory AND the generator can recover all three from a summary file name:
//   summary : etc.mu2e.VDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt
//   model   : nts.mu2e.VDResamplerModel_VD<id>_pdg<pdgTok>_<role>.<ver>.<run6>_<src8>.dat
//             (role: stage1/stage2/allAtOnce)
//   ROOT    : nts.mu2e.VDResamplerConfigure_VD<id>_hitDump.<ver>.<run6>_<src8>.root
//             (recommended TFileService name; module-external)
// <pdgTok> is 'm<abs>' for negative pdgIds; <id> is the VD id; <src8> is dataSourceIndex() of the
// dataSourceTag. The summary keeps its CSV *content* (comma-separated) under the .txt extension;
// the model keeps its binary content under .dat.
//
// versionTag MUST NOT contain '.' — the parser splits the name on dots, so a dotted version tag
// would be ambiguous. dataSourceTag MUST be one of the enumerated sources (see kDataSourceNames);
// an unrecognized tag is a hard error at encode time rather than a name the generator cannot parse.
//
// The producer (VDResamplerConfigure) and the consumer (VDResamplerGenerateMix) BOTH call these
// helpers, so a change here can never desynchronize the two sides. VDResamplerGenerateFromModel
// takes explicit model-file paths, so it does not use these builders.
//
// Header-only, dependency-free (just the standard library, plus cetlib_except for the encode-time
// errors). Yongyi Wu, Jul. 2026

#include "cetlib_except/exception.h"

#include <cstddef>
#include <string>
#include <vector>

namespace mu2e {
namespace VDResampler {

// Replace every character that is not [A-Za-z0-9_] with '_', so a tag is safe in a file name.
inline std::string sanitizeTag(const std::string& tag) {
  std::string s = tag;
  for (char& c : s) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
    if (!ok) c = '_';
  }
  return s;
}

// 'm<abs>' for negative pdgIds (bare number otherwise), matching the model-file convention.
inline std::string pdgFileToken(int pdgId) {
  return (pdgId < 0) ? "m" + std::to_string(-pdgId) : std::to_string(pdgId);
}

// The enumerated data sources, in the order that defines their numeric encoding. The INDEX of a
// tag in this list is what lands in the file name, so entries may be APPENDED but never reordered
// or removed — doing so would silently re-point every existing name at a different source.
inline const std::vector<std::string>& dataSourceNames() {
  static const std::vector<std::string> kNames = { "EleBeam", "MuBeam", "TargetStops1809", "Neutrals" };
  return kNames;
}

// Index of dataSourceTag in the enumeration, or -1 if it is not an enumerated source.
inline int dataSourceIndex(const std::string& dataSourceTag) {
  const std::vector<std::string>& names = dataSourceNames();
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == dataSourceTag) return static_cast<int>(i);
  return -1;
}

// Inverse of dataSourceIndex; empty string if the index is out of range.
inline std::string dataSourceFromIndex(int index) {
  const std::vector<std::string>& names = dataSourceNames();
  if (index < 0 || static_cast<std::size_t>(index) >= names.size()) return std::string();
  return names[index];
}

// Zero-pad a non-negative integer to `width` digits. Values that already exceed the width are
// emitted in full rather than truncated — a wrong-but-complete number beats a silently mangled one.
inline std::string zeroPad(int value, std::size_t width) {
  std::string s = std::to_string(value);
  return (s.size() >= width) ? s : std::string(width - s.size(), '0') + s;
}

// The <run6>_<src8> sequencer field shared by every name. Throws if dataSourceTag is not one of the
// enumerated sources, or if runNumber is negative.
inline std::string sequencerField(int runNumber, const std::string& dataSourceTag) {
  if (runNumber < 0)
    throw cet::exception("VDResamplerNameHelper")
      << "runNumber must be non-negative (got " << runNumber << ").";
  const int srcIndex = dataSourceIndex(dataSourceTag);
  if (srcIndex < 0) {
    cet::exception e("VDResamplerNameHelper");
    e << "dataSourceTag \"" << dataSourceTag << "\" is not an enumerated data source. Expected one of:";
    for (const std::string& n : dataSourceNames()) e << " " << n;
    e << ". Add it to dataSourceNames() (APPEND only) if it is a new source.";
    throw e;
  }
  return zeroPad(runNumber, 6) + "_" + zeroPad(srcIndex, 8);
}

// versionTag is embedded between two dots, so a dotted tag would break parseSummaryFileName.
inline std::string checkedVersionTag(const std::string& versionTag) {
  if (versionTag.find('.') != std::string::npos)
    throw cet::exception("VDResamplerNameHelper")
      << "versionTag \"" << versionTag << "\" must not contain '.'.";
  return sanitizeTag(versionTag);
}

inline std::string summaryFileName(const std::string& versionTag, int virtualDetectorID,
                                   const std::string& dataSourceTag, int runNumber) {
  return "etc.mu2e.VDResamplerConfigure_VD" + std::to_string(virtualDetectorID) + "_hitSummary."
       + checkedVersionTag(versionTag) + "." + sequencerField(runNumber, dataSourceTag) + ".txt";
}

// role is "stage1" | "stage2" | "allAtOnce". Extension is .dat (binary, full double precision).
inline std::string modelFileName(const std::string& role, const std::string& versionTag,
                                 int virtualDetectorID, const std::string& dataSourceTag,
                                 int pdgId, int runNumber) {
  return "nts.mu2e.VDResamplerModel_VD" + std::to_string(virtualDetectorID)
       + "_pdg" + pdgFileToken(pdgId) + "_" + role + "."
       + checkedVersionTag(versionTag) + "." + sequencerField(runNumber, dataSourceTag) + ".dat";
}

inline std::string rootDumpFileName(const std::string& versionTag, int virtualDetectorID,
                                    const std::string& dataSourceTag, int runNumber) {
  return "nts.mu2e.VDResamplerConfigure_VD" + std::to_string(virtualDetectorID) + "_hitDump."
       + checkedVersionTag(versionTag) + "." + sequencerField(runNumber, dataSourceTag) + ".root";
}

// Recover (versionTag, virtualDetectorID, dataSourceTag, runNumber) from a summary file name of the
// form ".../etc.mu2e.VDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt". Strips any
// directory prefix and returns false if the pattern does not match.
//
// The name is split on '.' into exactly five fields: tier, "mu2e", description, version, sequencer
// (the ".txt" extension is stripped first). That is unambiguous because versionTag is dot-free (see
// checkedVersionTag) — but the description field may itself contain underscores, so within it we
// anchor on the fixed "VDResamplerConfigure_VD" prefix and "_hitSummary" suffix. The source is
// recovered by decoding the 8-digit index back through the enumeration, so an index written by a
// newer build with extra sources will fail here rather than resolve to the wrong source.
inline bool parseSummaryFileName(const std::string& path, std::string& versionTag,
                                 int& virtualDetectorID, std::string& dataSourceTag,
                                 int& runNumber) {
  const std::size_t slash = path.find_last_of("/\\");
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);

  const std::string ext = ".txt";
  if (base.size() <= ext.size()) return false;
  if (base.compare(base.size() - ext.size(), ext.size(), ext) != 0) return false;
  base = base.substr(0, base.size() - ext.size());

  // Split the remainder on '.' — expect exactly {"etc", "mu2e", <desc>, <ver>, <seq>}.
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t dot = base.find('.', start);
    if (dot == std::string::npos) { fields.push_back(base.substr(start)); break; }
    fields.push_back(base.substr(start, dot - start));
    start = dot + 1;
  }
  if (fields.size() != 5) return false;
  if (fields[0] != "etc" || fields[1] != "mu2e") return false;

  // <desc> = "VDResamplerConfigure_VD<id>_hitSummary"
  const std::string pre = "VDResamplerConfigure_VD";
  const std::string suf = "_hitSummary";
  const std::string& desc = fields[2];
  if (desc.size() <= pre.size() + suf.size()) return false;
  if (desc.compare(0, pre.size(), pre) != 0) return false;
  if (desc.compare(desc.size() - suf.size(), suf.size(), suf) != 0) return false;
  const std::string idStr = desc.substr(pre.size(), desc.size() - pre.size() - suf.size());
  if (idStr.empty() || idStr.find_first_not_of("0123456789") != std::string::npos) return false;

  const std::string& ver = fields[3];
  if (ver.empty()) return false;

  // <seq> = "<run6>_<src8>"
  const std::string& seq = fields[4];
  const std::size_t us = seq.find('_');
  if (us == std::string::npos) return false;
  const std::string runStr = seq.substr(0, us);
  const std::string srcStr = seq.substr(us + 1);
  if (runStr.empty() || runStr.find_first_not_of("0123456789") != std::string::npos) return false;
  if (srcStr.empty() || srcStr.find_first_not_of("0123456789") != std::string::npos) return false;
  const std::string src = dataSourceFromIndex(std::stoi(srcStr));
  if (src.empty()) return false;

  versionTag        = ver;
  virtualDetectorID = std::stoi(idStr);
  dataSourceTag     = src;
  runNumber         = std::stoi(runStr);
  return true;
}

} // namespace VDResampler
} // namespace mu2e
