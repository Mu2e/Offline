// Counters and one end-of-job report for a DQM client.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMDiagnostics.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <sstream>

namespace mu2e {

void DQMDiagnostics::Count(const std::string& key, const std::string& message)
{
  const bool first = counts_.find(key) == counts_.end();
  ++counts_[key];
  if (first) {
    messages_[key] = message;
    mf::LogWarning(name_) << message << " Reported once per job; the total is in "
                          << "the end-of-job summary.";
  }
}

void DQMDiagnostics::Note(const std::string& key, double value)
{
  auto it = maxima_.find(key);
  if (it == maxima_.end() || value > it->second) {
    maxima_[key] = value;
  }
}

long long DQMDiagnostics::count(const std::string& key) const
{
  auto it = counts_.find(key);
  return it == counts_.end() ? 0 : it->second;
}

double DQMDiagnostics::max(const std::string& key) const
{
  auto it = maxima_.find(key);
  return it == maxima_.end() ? 0. : it->second;
}

void DQMDiagnostics::Report() const
{
  if (empty()) {
    return;
  }
  std::ostringstream out;
  out << "diagnostics summary:";
  for (const auto& [key, n] : counts_) {
    out << "\n  " << key << ": " << n;
    auto m = maxima_.find(key);
    if (m != maxima_.end()) {
      out << " (max " << m->second << ")";
    }
    auto message = messages_.find(key);
    if (message != messages_.end() && message->second != key) {
      out << " -- " << message->second;
    }
  }
  for (const auto& [key, value] : maxima_) {
    if (counts_.find(key) == counts_.end()) {
      out << "\n  " << key << ": max " << value;
    }
  }
  mf::LogWarning(name_) << out.str();
}

} // namespace mu2e
