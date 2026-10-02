#ifndef DQMHelpers_inc_DQMDiagnostics_hh
#define DQMHelpers_inc_DQMDiagnostics_hh
// Counters for the conditions a DQM client meets but does not histogram:
// an out-of-range id, a value past the end of a fixed axis, a missing partner.
//
// With fixed binning these are how an off-scale value gets reported at all, so
// they must reach the log rather than sit in an unread member. Count() warns
// the first time and totals are printed at end of job.
//
// Original Author: R. Mina

#include <map>
#include <string>

namespace mu2e {

class DQMDiagnostics {
public:
  void SetName(const std::string& name) { name_ = name; }

  // One occurrence. The message is logged once per key.
  void Count(const std::string& key, const std::string& message);
  void Count(const std::string& key) { Count(key, key); }
  // Largest value seen for a key, for "the axis ends at X but the data reached Y".
  void Note(const std::string& key, double value);

  long long count(const std::string& key) const;
  double max(const std::string& key) const;
  bool empty() const { return counts_.empty() && maxima_.empty(); }
  // One LogWarning listing every key with its count and maximum. Call at end
  // of job; does nothing when nothing was counted.
  void Report() const;

private:
  std::string name_{"DQMDiagnostics"};
  std::map<std::string, long long> counts_;
  std::map<std::string, std::string> messages_;
  std::map<std::string, double> maxima_;
};

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMDiagnostics_hh */
