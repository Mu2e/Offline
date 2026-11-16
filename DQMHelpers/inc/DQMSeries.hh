#ifndef DQMHelpers_inc_DQMSeries_hh
#define DQMHelpers_inc_DQMSeries_hh
// A capped time series (TGraph) for the online monitor.
//
// Graphs do not merge: hadd concatenates their points instead of adding a
// series, so a series is never part of the offline, mergeable output. A client
// books them only when its set's liveSeries flag is on, which offline presets
// leave off; with it off, book() hands back an inert series and every Add is a
// no-op, so the client fills unconditionally.
//
// Original Author: R. Mina

#include "art_root_io/TFileDirectory.h"

#include "TGraph.h"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace mu2e {

class DQMSeries {
public:
  DQMSeries() = default;
  DQMSeries(std::unique_ptr<TGraph> graph, std::size_t maxPoints);

  bool valid() const { return graph_ != nullptr; }
  TGraph* graph() const { return graph_.get(); }

  // Append a point, dropping the oldest once the cap is reached.
  void Add(double x, double y);
  // Append only when y changed.
  void AddIfChanged(double x, double y);
  // Step graph of a slowly changing value (a latency, a status word): a change
  // adds a vertical edge at x, an unchanged value extends the last step to x.
  void Step(double x, double y);
  void Clear();

private:
  void trim();

  std::unique_ptr<TGraph> graph_;
  std::size_t maxPoints_{0};
  std::optional<double> last_{};
};

// Books and owns the series of one client. Handles stay valid for its lifetime.
class DQMSeriesSet {
public:
  static constexpr std::size_t kDefaultMaxPoints = 10000;

  void Book(art::TFileDirectory dir, bool enabled);
  bool enabled() const { return enabled_; }

  // `path` may name a subdirectory, e.g. "graphs/g_linkLatency_link0".
  // maxPoints 0 takes kDefaultMaxPoints. Booking a path twice returns the first.
  DQMSeries& book(const std::string& path, const std::string& title,
                  std::size_t maxPoints = 0);
  DQMSeries& get(const std::string& path);
  // In booking order: what the online monitor ships.
  const std::vector<TGraph*>& graphs() const { return graphs_; }
  void ResetContents();
  // Write a copy of every graph at its path. Once per job.
  void Persist();

private:
  art::TFileDirectory& directoryFor(const std::string& dirPath);

  std::optional<art::TFileDirectory> dir_;
  std::map<std::string, art::TFileDirectory> subdirs_;
  bool enabled_{false};
  bool persisted_{false};
  std::map<std::string, std::unique_ptr<DQMSeries>> series_;
  std::vector<std::string> order_;
  std::vector<TGraph*> graphs_;
  DQMSeries inert_;
};

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMSeries_hh */
