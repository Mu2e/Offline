// FHiCL parsing for DQMSegmentation::Config.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMSegmentationConfig.hh"

#include "cetlib_except/exception.h"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace mu2e {

namespace {

void rejectUnknown(const fhicl::ParameterSet& ps,
                   const std::set<std::string>& allowed, const char* where)
{
  for (const auto& key : ps.get_names()) {
    if (allowed.find(key) == allowed.end()) {
      std::string known;
      for (const auto& a : allowed) {
        known += known.empty() ? a : ", " + a;
      }
      throw cet::exception("DQMSegmentation")
          << "unknown key \"" << key << "\" in " << where << ". Known keys: "
          << known << ".\n";
    }
  }
}

DQMSegmentation::SubRunConfig parseSubRun(const fhicl::ParameterSet& ps)
{
  rejectUnknown(ps, {"keep", "persist", "persistLive"}, "segmentation rule \"subrun\"");
  DQMSegmentation::SubRunConfig c;
  c.keep = ps.get<int>("keep", 0);
  c.persist = ps.get<bool>("persist", true);
  c.persistLive = ps.get<bool>("persistLive", false);
  return c;
}

DQMSegmentation::WindowConfig parseWindow(const fhicl::ParameterSet& ps)
{
  rejectUnknown(ps,
                {"span", "unit", "subdivisions", "keep", "persist", "persistLive"},
                "segmentation rule \"window\"");
  DQMSegmentation::WindowConfig c;
  c.spanSet = ps.has_key("span");
  const long long span = ps.get<long long>("span", 50000);
  if (span < 1) {
    throw cet::exception("DQMSegmentation")
        << "segmentation window span must be >= 1, got " << span << ".\n";
  }
  c.span = static_cast<unsigned long long>(span);

  const std::string unit = ps.get<std::string>("unit", "event");
  if (!DQMSegmentation::unitFromString(unit, c.unit)) {
    throw cet::exception("DQMSegmentation")
        << "unknown segmentation window unit \"" << unit
        << "\". Known units: event, ewt, subrun.\n";
  }

  c.subdivisions = std::max(ps.get<int>("subdivisions", 10), 1);
  c.keep = std::max(ps.get<int>("keep", 0), 0);
  c.persist = ps.get<bool>("persist", false);
  c.persistLive = ps.get<bool>("persistLive", false);
  return c;
}

DQMSegmentation::Rule parseRule(const fhicl::ParameterSet& ps)
{
  rejectUnknown(ps,
                {"match", "enabled", "modes", "job", "subrun", "window", "liveName",
                 "publish", "group"},
                "a segmentation rule");

  DQMSegmentation::Rule rule;
  rule.match = ps.get<std::string>("match", "*");
  rule.enabled = ps.get<bool>("enabled", true);
  rule.liveName = ps.get<std::string>("liveName", "");
  rule.publish = ps.get<bool>("publish", false);
  rule.group = ps.get<std::string>("group", "");

  const auto modes = ps.get<std::vector<std::string>>(
      "modes", std::vector<std::string>{"job"});
  rule.job = false;
  for (const auto& mode : modes) {
    if (mode == "job") {
      rule.job = true;
    } else if (mode == "subrun") {
      rule.subrun.enabled = true;
    } else if (mode == "window") {
      rule.window.enabled = true;
    } else {
      throw cet::exception("DQMSegmentation")
          << "unknown segmentation mode \"" << mode
          << "\" for match \"" << rule.match
          << "\". Known modes: job, subrun, window.\n";
    }
  }

  fhicl::ParameterSet sub;
  if (ps.get_if_present<fhicl::ParameterSet>("job", sub)) {
    rejectUnknown(sub, {"persist"}, "segmentation rule \"job\"");
    rule.jobPersist = sub.get<bool>("persist", true);
  }
  if (ps.get_if_present<fhicl::ParameterSet>("subrun", sub)) {
    const bool enabled = rule.subrun.enabled;
    rule.subrun = parseSubRun(sub);
    rule.subrun.enabled = enabled;
  }
  if (ps.get_if_present<fhicl::ParameterSet>("window", sub)) {
    const bool enabled = rule.window.enabled;
    rule.window = parseWindow(sub);
    rule.window.enabled = enabled;
  }

  if (!rule.job && !rule.subrun.enabled && !rule.window.enabled && rule.enabled) {
    throw cet::exception("DQMSegmentation")
        << "segmentation rule for match \"" << rule.match
        << "\" selects no modes. Use enabled: false to drop the histogram "
        << "instead.\n";
  }
  return rule;
}

} // namespace

DQMSegmentation::Config parseSegmentation(const fhicl::ParameterSet& ps)
{
  rejectUnknown(
      ps, {"annotateTitles", "stampMetadata", "subrunDir", "segmentDir", "rules"},
      "the segmentation block");

  DQMSegmentation::Config config;
  config.annotateTitles = ps.get<bool>("annotateTitles", true);
  config.stampMetadata = ps.get<bool>("stampMetadata", true);
  config.subrunDir = ps.get<std::string>("subrunDir", "bySubrun");
  config.segmentDir = ps.get<std::string>("segmentDir", "segments");

  const auto rules = ps.get<std::vector<fhicl::ParameterSet>>(
      "rules", std::vector<fhicl::ParameterSet>{});
  config.rules.reserve(rules.size());
  for (const auto& rule : rules) {
    config.rules.push_back(parseRule(rule));
  }
  return config;
}

DQMSegmentation::Config parseSegmentation(const fhicl::OptionalDelegatedParameter& p)
{
  fhicl::ParameterSet ps;
  if (!p.get_if_present<fhicl::ParameterSet>(ps)) {
    return DQMSegmentation::Config{};
  }
  return parseSegmentation(ps);
}

} // namespace mu2e
