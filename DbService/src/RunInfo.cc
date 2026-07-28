#include "Offline/DbService/inc/RunInfo.hh"

#include "Offline/GeneralUtilities/inc/TimeUtility.hh"

#include "cetlib_except/exception.h"
#include "nlohmann/json.hpp"

#include <ctime>
#include <string>

using namespace mu2e;

//**************************************************
// Collect the cat-2 DbService CID table entries for this run.  Only the "TRG"
// config holds the CID tables, so restrict the lookup to that config.  This
// avoids emitting multiple concatenated (and hence invalid) JSON documents.
std::string RunInfo::dbTables2(bool qjson) const {
  for (auto const& config : _configs) {
    if (config.subsystem() != "TRG") {
      continue;
    }
    std::string tables = config.dbTables2(qjson);
    if (!tables.empty()) {
      return tables;
    }
  }

  // No TRG config, or it had no CID entries.  Return an empty JSON array when
  // JSON output was requested, otherwise an empty string.
  return qjson ? "[]" : std::string();
}

//**************************************************
// Collect the cat-3 DbService table key-value pairs across every config for
// this run.  Several configs can contribute pairs, so when JSON output is
// requested the individual per-config objects are merged into a single JSON
// object (an empty object "{}" if nothing is found).
std::string RunInfo::dbTables3(bool qjson) const {
  if (qjson) {
    nlohmann::json merged = nlohmann::json::object();
    for (auto const& config : _configs) {
      std::string tables = config.dbTables3(qjson);
      if (tables.empty()) {
        continue;
      }
      nlohmann::json part = nlohmann::json::parse(tables, nullptr, false);
      if (part.is_discarded() || !part.is_object()) {
        continue;
      }
      for (auto it = part.begin(); it != part.end(); ++it) {
        if (merged.contains(it.key())) {
          throw cet::exception("RUNINFO_DUPLICATE_TABLE")
              << "RunInfo::dbTables3() found duplicate table key \""
              << it.key() << "\" across run configs\n";
        }
        merged[it.key()] = it.value();
      }
    }
    return merged.dump(1, ' ');
  }

  // Non-JSON: concatenate the per-config values, one per line.
  std::string result;
  for (auto const& config : _configs) {
    std::string tables = config.dbTables3(qjson);
    if (!tables.empty()) {
      result += tables;
      result += "\n";
    }
  }
  if (!result.empty() && result.back() == '\n') {
    result.pop_back();
  }
  return result;
}

//**************************************************
int RunInfo::isDone(int timeout_seconds) const {
  // Check if we have the required data
  if (_transitions.empty() || _subruns.empty()) {
    throw cet::exception("RUNINFO_INCOMPLETE")
        << "RunInfo::isDone() cannot determine status: "
        << "transitions or subruns list is empty\n";
  }

  // Get the last transition
  const RunTransition& lastTransition = _transitions.back();
  int lastTransitionType = lastTransition.typeId();

  // Check for user stop (transitions 0, 1, 6, 7)
  if (lastTransitionType == 0 || lastTransitionType == 1 ||
      lastTransitionType == 6 || lastTransitionType == 7) {
    return 1;  // user stop
  }

  // Check for error stop (transition 2)
  if (lastTransitionType == 2) {
    return 2;  // error stop
  }

  // Check for timeout (probable crash)
  // Get the last transition time
  std::string lastTransTime = lastTransition.transitionTime();
  std::time_t lastTransTimeT;
  int rc = TimeUtility::parseTimeTZ(lastTransTime, lastTransTimeT);
  if (rc != 0) {
    throw cet::exception("RUNINFO_BAD_TIME")
        << "RunInfo::isDone() cannot parse transition time: " << lastTransTime
        << "\n";
  }

  // Get the last subrun stop time
  const RunSubRun& lastSubrun = _subruns.back();
  std::time_t lastSubrunTimeT = lastSubrun.stopTimeUnix();

  // Use the more recent of the two times
  std::time_t lastActivityTime = std::max(lastTransTimeT, lastSubrunTimeT);

  // Get current time
  std::time_t currentTime = std::time(nullptr);

  // Check if timeout exceeded
  if ((currentTime - lastActivityTime) > timeout_seconds) {
    return 3;  // probably crash end
  }

  // Still running
  return 0;
}
