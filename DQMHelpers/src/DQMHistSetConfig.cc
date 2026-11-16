// FHiCL parsing for DQMHistSet::Config.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMHistSetConfig.hh"

#include "cetlib_except/exception.h"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace mu2e {

DQMHistSet::Config toConfig(const fhicl::ParameterSet& ps)
{
  return toConfig(fhicl::Table<DQMHistSetFhicl>(ps)());
}

DQMHistSet::Config toConfig(const DQMHistSetFhicl& c)
{
  DQMHistSet::Config out;
  out.annotateTitles = c.annotateTitles();
  out.stampMetadata = c.stampMetadata();
  out.subrunDir = c.subrunDir();
  out.segmentDir = c.segmentDir();
  out.liveSeries = c.liveSeries();

  //Table<S>::value_type is S, so a sequence of tables reads back as the structs
  std::vector<DQMRuleFhicl> rules;
  if (!c.rules(rules)) {
    return out;  //no rules: job-only, the default
  }

  for (const DQMRuleFhicl& r : rules) {
    DQMHistSet::Rule rule;
    rule.match = r.match();
    rule.liveName = r.liveName();
    rule.publish = r.publish();
    rule.group = r.group();
    rule.archiveGroup = r.archiveGroup();
    rule.jobPersist = r.jobPersist();

    rule.job = false;
    for (const std::string& mode : r.modes()) {
      if (mode == "job") {
        rule.job = true;
      } else if (mode == "subrun") {
        rule.subrun.enabled = true;
      } else if (mode == "window") {
        rule.window.enabled = true;
      } else {
        throw cet::exception("DQMHistSetConfig")
            << "unknown mode \"" << mode << "\" in rule \"" << rule.match
            << "\"; expected job, subrun or window.\n";
      }
    }
    if (!rule.job && !rule.subrun.enabled && !rule.window.enabled) {
      throw cet::exception("DQMHistSetConfig")
          << "rule \"" << rule.match << "\" selects no modes.\n";
    }

    DQMSubRunFhicl subrun;
    if (r.subrun(subrun)) {
      rule.subrun.keep = subrun.keep();
      rule.subrun.persist = subrun.persist();
      rule.subrun.persistLive = subrun.persistLive();
    }

    DQMWindowFhicl window;
    if (r.window(window)) {
      rule.window.span = window.span();
      rule.window.spanSet = true;
      rule.window.subdivisions = std::max(window.subdivisions(), 1);
      rule.window.keep = window.keep();
      rule.window.persist = window.persist();
      rule.window.persistLive = window.persistLive();
      if (!DQMHistSet::unitFromString(window.unit(), rule.window.unit)) {
        throw cet::exception("DQMHistSetConfig")
            << "unknown window unit \"" << window.unit() << "\" in rule \""
            << rule.match << "\"; expected event, ewt or subrun.\n";
      }
    }

    out.rules.push_back(rule);
  }
  return out;
}

} // namespace mu2e
