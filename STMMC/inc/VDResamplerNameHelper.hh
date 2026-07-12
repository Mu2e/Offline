#pragma once

// File-name convention for the VD resampler pipeline (single source of truth).
//
// Every generated artifact embeds versionTag and dataSourceTag so several campaigns can coexist
// in one directory AND the generator can recover both from a summary file name:
//   summary : VDResamplerConfigure_<ver>_VD<id>_<src>_HitSummary.csv
//   model   : SBDMmodel_<role>_<ver>_VD<id>_<src>_pdg<pdgTok>.bin   (role: stage1/stage2/allAtOnce)
//   ROOT    : VDResamplerDump_<ver>_VD<id>_<src>.root   (recommended TFileService name; module-external)
// <src> is the sanitized dataSourceTag; <pdgTok> is 'm<abs>' for negative pdgIds; <id> is the VD id.
//
// The producer (VDResamplerConfigure) and the consumer (VDResamplerGenerateMix) BOTH call these
// helpers, so a change here can never desynchronize the two sides. VDResamplerGenerateFromModel
// takes explicit model-file paths, so it does not use these builders.
//
// Header-only, dependency-free (just the standard library). Yongyi Wu, Jul. 2026

#include <cstddef>
#include <string>

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

inline std::string summaryFileName(const std::string& versionTag, int virtualDetectorID,
                                   const std::string& dataSourceTag) {
  return "VDResamplerConfigure_" + sanitizeTag(versionTag) + "_VD" + std::to_string(virtualDetectorID)
       + "_" + sanitizeTag(dataSourceTag) + "_HitSummary.csv";
}

// role is "stage1" | "stage2" | "allAtOnce". Extension is .bin (full double precision).
inline std::string modelFileName(const std::string& role, const std::string& versionTag,
                                 int virtualDetectorID, const std::string& dataSourceTag, int pdgId) {
  return "SBDMmodel_" + role + "_" + sanitizeTag(versionTag) + "_VD" + std::to_string(virtualDetectorID)
       + "_" + sanitizeTag(dataSourceTag) + "_pdg" + pdgFileToken(pdgId) + ".bin";
}

inline std::string rootDumpFileName(const std::string& versionTag, int virtualDetectorID,
                                    const std::string& dataSourceTag) {
  return "VDResamplerDump_" + sanitizeTag(versionTag) + "_VD" + std::to_string(virtualDetectorID)
       + "_" + sanitizeTag(dataSourceTag) + ".root";
}

// Recover (versionTag, virtualDetectorID, dataSourceTag) from a summary file name of the form
// ".../VDResamplerConfigure_<ver>_VD<id>_<src>_HitSummary.csv". Strips any directory prefix and
// returns false if the pattern does not match. We anchor on the fixed "VDResamplerConfigure_"
// prefix and "_HitSummary.csv" suffix; the version is the FIRST underscore-delimited field, then a
// "VD<id>" field, then the remainder (which may contain underscores) is the source. versionTag is
// not expected to contain underscores; if it does, this split is ambiguous — keep versionTag simple.
inline bool parseSummaryFileName(const std::string& path, std::string& versionTag,
                                 int& virtualDetectorID, std::string& dataSourceTag) {
  const std::size_t slash = path.find_last_of("/\\");
  const std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  const std::string pre = "VDResamplerConfigure_";
  const std::string suf = "_HitSummary.csv";
  if (base.size() <= pre.size() + suf.size()) return false;
  if (base.compare(0, pre.size(), pre) != 0) return false;
  if (base.compare(base.size() - suf.size(), suf.size(), suf) != 0) return false;
  const std::string mid = base.substr(pre.size(), base.size() - pre.size() - suf.size()); // "<ver>_VD<id>_<src>"

  const std::size_t us1 = mid.find('_');                       // end of <ver>
  if (us1 == std::string::npos || us1 == 0) return false;
  const std::size_t us2 = mid.find('_', us1 + 1);              // end of VD<id>
  if (us2 == std::string::npos || us2 + 1 >= mid.size()) return false;
  const std::string ver = mid.substr(0, us1);
  const std::string vd  = mid.substr(us1 + 1, us2 - us1 - 1);  // "VD<id>"
  const std::string src = mid.substr(us2 + 1);
  if (vd.size() <= 2 || vd.compare(0, 2, "VD") != 0) return false;
  const std::string idStr = vd.substr(2);
  if (idStr.empty() || idStr.find_first_not_of("0123456789") != std::string::npos) return false;
  versionTag        = ver;
  virtualDetectorID = std::stoi(idStr);
  dataSourceTag     = src;
  return true;
}

} // namespace VDResampler
} // namespace mu2e
