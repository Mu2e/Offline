#ifndef DQMHelpers_inc_DQMHistSet_hh
#define DQMHelpers_inc_DQMHistSet_hh
// FHiCL-driven histogram segmentation for the shared DQM helpers.
//
// A helper books every histogram through this registry instead of calling
// art::TFileDirectory::make directly. The registry consults a rule list -- set
// from FHiCL, so it can change between runs with no rebuild -- and creates, per
// histogram, any of:
//
//   job     the histogram as it has always been: one copy, never reset
//   subrun  a live copy for the subrun in progress, plus archived copies of
//           previous subruns, named by real run/subrun so they survive hadd
//   window  a rolling copy covering the last N events / EWTs / subruns, plus
//           archived copies of previous spans
//
// The default is job-only, which reproduces the pre-segmentation output exactly.
//
// The window is a ring of `subdivisions` sub-blocks of span/subdivisions units
// each. Entries go into the current sub-block and into the live copy, so the
// live copy is exact and never lags; on rotation the live copy is rebuilt by
// summing the ring, which keeps bin errors and the entry count right. It
// therefore covers between span*(S-1)/S and span units -- the realised range is
// stamped in every copy's title, so no copy's meaning is ambiguous.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMAxis.hh"
#include "Offline/DQMHelpers/inc/DQMHist.hh"

#include "art_root_io/TFileDirectory.h"

#include "TH1.h"
#include "TH1F.h"

#include <cstdint>
#include <functional>
#include <ostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace mu2e {

class DQMHistSet {
public:
  enum class Unit { Event, Ewt, SubRun };

  static const char* unitName(Unit u);
  static bool unitFromString(const std::string& s, Unit& u);

  struct SubRunConfig {
    bool enabled{false};
    int keep{0};              //archived subruns to retain; -1 = all
    bool persist{true};       //write the archived copies
    bool persistLive{false};  //write the live (in-progress) copy
  };

  struct WindowConfig {
    bool enabled{false};
    unsigned long long span{50000};
    bool spanSet{false};      //false: the caller's own default may override
    Unit unit{Unit::Event};
    int subdivisions{10};     //ring depth; 1 = disjoint blocks
    int keep{0};              //archived previous spans to retain
    bool persist{false};      //write the archived copies
    bool persistLive{false};  //write the live copy
  };

  // There is no "enabled": FHiCL chooses copies of a histogram, never whether
  // it exists, so every job of a client writes the same set.
  struct Rule {
    std::string match{"*"};   //glob over the directory-qualified histogram path
    bool job{true};
    bool jobPersist{true};
    SubRunConfig subrun{};
    WindowConfig window{};
    std::string liveName{};   //explicit name for the window live copy
    // Publishing. The registry does not know what a consumer does with these:
    // `group` is an opaque label, and it is the caller that decides whether it
    // means an otsdaq HistoSender folder, a web tab or anything else. Keeping
    // it opaque is what lets DQMHelpers build in the DAQ process without
    // knowing HistoSender exists.
    bool publish{false};      //include this histogram's copies in publishedCopies()
    std::string group{};      //label to collect them under; empty = each copy's own name
    std::string archiveGroup{};  //label for the _prevN and per-subrun archives; empty = `group`
  };

  struct Config {
    bool annotateTitles{true};
    bool stampMetadata{true};
    std::string subrunDir{"bySubrun"};
    std::string segmentDir{"segments"};
    // Refuse a book*() after FreezeBooking(), or a second booking of one path:
    // the histogram set must not depend on what the input happens to contain.
    // DQMClient turns both of these on.
    bool strictBooking{false};
    // Book and fill "nEvents" here rather than in each client.
    bool autoNEvents{false};
    // Book the client's DQMSeries (online only; graphs do not merge).
    bool liveSeries{false};
    std::vector<Rule> rules{};
  };

  // What one copy covers. Written into the title and the dqmSegment stamp.
  struct Range {
    std::size_t firstEvent{0};
    std::size_t lastEvent{0};
    std::size_t nEvents{0};
    unsigned long long firstClock{0};
    unsigned long long lastClock{0};
    int run{-1};
    int subrun{-1};
    bool complete{false};
    bool started{false};
  };

  DQMHistSet() = default;
  explicit DQMHistSet(const Config& config) : config_(config) {}

  void SetConfig(const Config& config) { config_ = config; }
  const Config& config() const { return config_; }

  // Call once before any book*(). Everything is created under `dir`.
  void Book(art::TFileDirectory dir);
  bool booked() const { return dir_.has_value(); }

  // Identity of the client filling this set, written into the output as
  // dqmBinningVersion. A merging or metrics tool refuses to combine versions.
  void SetVersion(const std::string& clientName, int binningVersion);
  // No histogram may be booked after this; see Config::strictBooking.
  void FreezeBooking() { frozen_ = true; }
  // name, type and axes of every booked histogram, one per line: the reference
  // that a fixed-binning check compares against.
  void WriteCatalogue(std::ostream& out) const;

  // `path` may name a subdirectory, e.g. "timing/dtFpgaPairs". Rules match
  // against that whole path. Binning comes only from DQMAxis constants, so it
  // is reviewable in the client header and dumped by WriteCatalogue.
  template <class H>
  DQMH1<H> book1(const std::string& path, const std::string& title,
                 const DQMAxis& x);
  template <class H>
  DQMH2<H> book2(const std::string& path, const std::string& title,
                 const DQMAxis& x, const DQMAxis& y);
  // A summary filled once at end of job from other histograms (a fit result, a
  // rate). It takes no rule: segment copies of it would all be identical.
  template <class H>
  DQMH1<H> bookSummary1(const std::string& path, const std::string& title,
                        const DQMAxis& x);
  template <class H>
  DQMH2<H> bookSummary2(const std::string& path, const std::string& title,
                        const DQMAxis& x, const DQMAxis& y);

  // Clock hooks. Advance() must be called once per event, before any fill, so
  // the window ring rotates between events rather than inside one.
  void BeginSubRun(int run, int subrun);
  void EndSubRun();
  void Advance(std::size_t eventIndex, std::optional<uint64_t> ewt);
  // Rewrite the titles and stamps of the live (in-progress) copies from the
  // range they currently hold. Rotations and Finalize() do this already; call
  // it directly before reading a live copy's title mid-span, as the online
  // monitor does before it ships or draws one.
  void RefreshLabels();
  void Finalize();
  // Empty every copy and restart the window ring, keeping the objects booked
  // (online, a new run must not inherit the last one). Named per-subrun
  // archives describe finished subruns, so they are kept.
  void ResetContents();

  // Every copy of one histogram, job copy first: what the online monitor ships
  // so that a newly configured segment needs no C++ change to reach the GUI.
  std::vector<TH1*> copies(const std::string& path) const;
  std::vector<TH1*> allCopies() const;
  // The rolling window copy, or nullptr when this histogram has no window mode.
  TH1* live(const std::string& path) const;

  // Every copy of every histogram whose rule set `publish`, collected under
  // that rule's `group`. A rule with no group gives each copy an entry of its
  // own, keyed on the copy's object name -- so a live window copy named by
  // `liveName` is published under that name, which is what a GUI subscribing
  // to a fixed name wants. This is the whole of what a consumer needs: which
  // histograms to ship and how to group them, both chosen in FHiCL.
  std::map<std::string, std::vector<TH1*>> publishedCopies() const;

  std::size_t nHistograms() const { return entries_.size(); }
  // Set once if an "ewt" rule ran on input with no event window tag.
  bool ewtFallback() const { return ewtFallback_; }

private:
  using DirFactory =
      std::function<TH1*(art::TFileDirectory&, const std::string&, const std::string&)>;
  using OwnedFactory =
      std::function<std::unique_ptr<TH1>(const std::string&, const std::string&)>;

  struct Entry {
    std::string path;
    std::string base;
    std::string dirPath;
    std::string title;
    Rule rule;

    DQMHistTargets targets;
    DirFactory makeInDir;
    OwnedFactory makeOwned;

    TH1* job{nullptr};
    TH1* subrunLive{nullptr};
    TH1* windowLive{nullptr};
    std::vector<TH1*> subBlocks;
    std::vector<TH1*> windowArchive;  //_prev1 .. _prevK, newest first
    std::vector<TH1*> subrunArchive;  //created at each subrun end

    Range jobRange;
    Range subrunRange;
    Range liveRange;
    Range spanRange;                  //the span the ring is accumulating
    std::vector<Range> subBlockRange;
    std::vector<Range> archiveRange;

    bool windowStarted{false};
    unsigned long long subBlockStartClock{0};
    unsigned long long subBlockWidth{1};
    std::size_t cur{0};
    std::size_t completedSubBlocks{0};
    std::size_t nSubRunsArchived{0};
  };

  static bool globMatch(const std::string& pattern, const std::string& text);
  static std::string annotate(const std::string& title, const std::string& tag);
  static std::string stampText(const Range& range, const char* mode, int index);

  const Rule& ruleFor(const std::string& path) const;
  art::TFileDirectory& directoryFor(const std::string& dirPath);
  TH1* create(Entry& entry, const std::string& dirPath, const std::string& name,
              bool persist);
  Entry* addEntry(const std::string& path, const std::string& title,
                  DirFactory makeInDir, OwnedFactory makeOwned, bool summary);
  template <class H, class... Args>
  Entry* bookEntry(const std::string& path, const std::string& title,
                   bool summary, Args... args);
  void refreshTargets(Entry& entry);
  void dropOwned(TH1* h);
  static void note(Range& range, std::size_t event, unsigned long long clock,
                   int run, int subrun);
  void advanceWindow(Entry& entry, unsigned long long clock);
  void rotateSubBlock(Entry& entry, unsigned long long clock);
  void archiveSpan(Entry& entry);
  void rebuildLive(Entry& entry);
  void label(Entry& entry, TH1* h, const Range& range, const char* mode,
             int index);
  void labelAll(Entry& entry);
  unsigned long long clockFor(const Entry& entry) const;

  Config config_;
  std::optional<art::TFileDirectory> dir_;
  std::map<std::string, art::TFileDirectory> subdirs_;
  std::vector<std::unique_ptr<Entry>> entries_;
  // Entries with a subrun or window copy. Advance() runs once per event and
  // there can be thousands of job-only per-FPGA histograms, none of which
  // need per-event bookkeeping.
  std::vector<Entry*> segmented_;
  std::map<std::string, Entry*> byPath_;
  std::vector<std::unique_ptr<TH1>> owned_;

  std::size_t currentEvent_{0};
  unsigned long long currentEwt_{0};
  bool haveEwt_{false};
  bool ewtFallback_{false};
  int run_{-1};
  int subrun_{-1};
  unsigned long long subrunClock_{0};
  Rule defaultRule_{};
  bool frozen_{false};
  std::string clientName_{};
  int binningVersion_{0};
  DQMH1<TH1F> h_nEvents_;
};

template <class H, class... Args>
DQMHistSet::Entry* DQMHistSet::bookEntry(const std::string& path,
                                         const std::string& title, bool summary,
                                         Args... args)
{
  auto inDir = [args...](art::TFileDirectory& d, const std::string& n,
                         const std::string& t) -> TH1* {
    return d.make<H>(n.c_str(), t.c_str(), args...);
  };
  auto owned = [args...](const std::string& n,
                         const std::string& t) -> std::unique_ptr<TH1> {
    auto h = std::make_unique<H>(n.c_str(), t.c_str(), args...);
    h->SetDirectory(nullptr);
    return h;
  };
  return addEntry(path, title, inDir, owned, summary);
}

template <class H>
DQMH1<H> DQMHistSet::book1(const std::string& path, const std::string& title,
                           const DQMAxis& x)
{
  Entry* e = bookEntry<H>(path, title, false, x.n, x.lo, x.hi);
  return e == nullptr ? DQMH1<H>() : DQMH1<H>(&e->targets, static_cast<H*>(e->job));
}

template <class H>
DQMH2<H> DQMHistSet::book2(const std::string& path, const std::string& title,
                           const DQMAxis& x, const DQMAxis& y)
{
  Entry* e = bookEntry<H>(path, title, false, x.n, x.lo, x.hi, y.n, y.lo, y.hi);
  return e == nullptr ? DQMH2<H>() : DQMH2<H>(&e->targets, static_cast<H*>(e->job));
}

template <class H>
DQMH1<H> DQMHistSet::bookSummary1(const std::string& path, const std::string& title,
                                  const DQMAxis& x)
{
  Entry* e = bookEntry<H>(path, title, true, x.n, x.lo, x.hi);
  return e == nullptr ? DQMH1<H>() : DQMH1<H>(&e->targets, static_cast<H*>(e->job));
}

template <class H>
DQMH2<H> DQMHistSet::bookSummary2(const std::string& path, const std::string& title,
                                  const DQMAxis& x, const DQMAxis& y)
{
  Entry* e = bookEntry<H>(path, title, true, x.n, x.lo, x.hi, y.n, y.lo, y.hi);
  return e == nullptr ? DQMH2<H>() : DQMH2<H>(&e->targets, static_cast<H*>(e->job));
}

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMHistSet_hh */
