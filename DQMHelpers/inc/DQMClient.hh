#ifndef DQMHelpers_inc_DQMClient_hh
#define DQMHelpers_inc_DQMClient_hh
// Base class for a subdetector's DQM client.
//
// A client owns a DQMHistSet, a DQMSeriesSet and a DQMDiagnostics, and the
// lifecycle below. A subdetector writes only book() and its own typed Fill(),
// which calls beginEvent() first -- the client is the only thing that knows
// where its clock comes from.
//
//   class TrkDigiDQM : public DQMClient {
//    public:
//     static constexpr int kBinningVersion = 1;
//     static constexpr DQMAxis kNDigis = DQMAxis::Counts(0, 100);
//     explicit TrkDigiDQM(const DQMHistSet::Config& hists = {}) :
//         DQMClient("TrkDigiDQM", kBinningVersion, hists) {}
//     void Fill(const StrawDigiCollection& digis) {
//       beginEvent();
//       nDigis_.Fill(digis.size());
//     }
//    private:
//     void book() override {
//       nDigis_ = hists().book1<TH1F>("nDigis", "Digis / event;N;Events", kNDigis);
//     }
//     DQMH1<TH1F> nDigis_;
//   };
//
// The module side is the same for every subdetector: Book() in beginJob,
// BeginSubRun / EndSubRun, Fill in analyze, EndJob in endJob, and
// ResetForNewRun in beginRun online. The set books and fills nEvents itself,
// and refuses a booking outside book() or a second booking of one name.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMDiagnostics.hh"
#include "Offline/DQMHelpers/inc/DQMHistSet.hh"
#include "Offline/DQMHelpers/inc/DQMSeries.hh"

#include "art_root_io/TFileDirectory.h"

#include <cstdint>
#include <optional>
#include <string>

namespace mu2e {

class DQMClient {
public:
  DQMClient(const std::string& name, int binningVersion,
            const DQMHistSet::Config& hists);
  virtual ~DQMClient() = default;

  const std::string& name() const { return name_; }
  int binningVersion() const { return binningVersion_; }

  // Books everything under `dir`, then freezes: a client's histogram set does
  // not depend on the input.
  void Book(art::TFileDirectory dir);
  bool booked() const { return booked_; }
  void BeginSubRun(int run, int subrun);
  void EndSubRun();
  // Runs the client's end-of-job work, finalizes the segment labels, writes
  // the series and prints the diagnostics summary.
  void EndJob();
  // Online: a new run must not inherit the last one's contents. nEvents()
  // keeps counting.
  void ResetForNewRun();

  DQMHistSet& hists() { return hists_; }
  const DQMHistSet& hists() const { return hists_; }
  DQMSeriesSet& series() { return series_; }
  const DQMSeriesSet& series() const { return series_; }
  DQMDiagnostics& diag() { return diag_; }
  const DQMDiagnostics& diag() const { return diag_; }
  std::size_t nEvents() const { return nEvents_; }
  int run() const { return run_; }
  int subrun() const { return subrun_; }

protected:
  // Book every histogram of this client. Called once, by Book().
  virtual void book() = 0;
  // Optional hooks, each called before the core's own step.
  virtual void endSubRun() {}
  virtual void endJob() {}
  virtual void resetForNewRun() {}
  // Call first in Fill(): counts the event and advances the window clock.
  void beginEvent(std::optional<uint64_t> clock = std::nullopt);

private:
  std::string name_;
  int binningVersion_{0};
  DQMHistSet hists_;
  DQMSeriesSet series_;
  DQMDiagnostics diag_;
  std::size_t nEvents_{0};
  int run_{-1};
  int subrun_{-1};
  bool booked_{false};
  bool ended_{false};
};

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMClient_hh */
