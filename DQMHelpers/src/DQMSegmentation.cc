// FHiCL-driven histogram segmentation for the shared DQM helpers.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMSegmentation.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TList.h"
#include "TNamed.h"
#include "TString.h"

#include <algorithm>
#include <sstream>

namespace mu2e {

const char* DQMSegmentation::unitName(Unit u)
{
  switch (u) {
    case Unit::Event:
      return "events";
    case Unit::Ewt:
      return "EWT";
    case Unit::SubRun:
      return "subruns";
  }
  return "events";
}

bool DQMSegmentation::unitFromString(const std::string& s, Unit& u)
{
  if (s == "event" || s == "events") {
    u = Unit::Event;
    return true;
  }
  if (s == "ewt" || s == "EWT") {
    u = Unit::Ewt;
    return true;
  }
  if (s == "subrun" || s == "subruns") {
    u = Unit::SubRun;
    return true;
  }
  return false;
}

// '*' matches any run of characters, '?' exactly one. Backtracking on the last
// '*' keeps it linear in practice for the patterns a rule list carries.
bool DQMSegmentation::globMatch(const std::string& pattern,
                                const std::string& text)
{
  std::size_t p = 0, t = 0, star = std::string::npos, mark = 0;
  while (t < text.size()) {
    if (p < pattern.size() && (pattern[p] == '?' || pattern[p] == text[t])) {
      ++p;
      ++t;
    } else if (p < pattern.size() && pattern[p] == '*') {
      star = p++;
      mark = t;
    } else if (star != std::string::npos) {
      p = star + 1;
      t = ++mark;
    } else {
      return false;
    }
  }
  while (p < pattern.size() && pattern[p] == '*') {
    ++p;
  }
  return p == pattern.size();
}

// ROOT titles carry the axis labels after the first ';', so the tag has to go
// in front of it or the x-axis label silently becomes part of the title.
std::string DQMSegmentation::annotate(const std::string& title,
                                      const std::string& tag)
{
  if (tag.empty()) {
    return title;
  }
  const std::size_t semi = title.find(';');
  if (semi == std::string::npos) {
    return title + " " + tag;
  }
  return title.substr(0, semi) + " " + tag + title.substr(semi);
}

std::string DQMSegmentation::stampText(const Range& range, const char* mode,
                                       int index)
{
  std::ostringstream os;
  os << "mode=" << mode << ";index=" << index << ";firstEvent="
     << range.firstEvent << ";lastEvent=" << range.lastEvent
     << ";nEvents=" << range.nEvents << ";firstClock=" << range.firstClock
     << ";lastClock=" << range.lastClock << ";run=" << range.run
     << ";subrun=" << range.subrun << ";complete=" << (range.complete ? 1 : 0);
  return os.str();
}

void DQMSegmentation::note(Range& range, std::size_t event,
                           unsigned long long clock, int run, int subrun)
{
  if (!range.started) {
    range.started = true;
    range.firstEvent = event;
    range.firstClock = clock;
    range.run = run;
    range.subrun = subrun;
  }
  range.lastEvent = event;
  range.lastClock = clock;
  ++range.nEvents;
}

void DQMSegmentation::Book(art::TFileDirectory dir)
{
  dir_ = dir;
}

art::TFileDirectory& DQMSegmentation::directoryFor(const std::string& dirPath)
{
  if (dirPath.empty()) {
    return *dir_;
  }
  auto it = subdirs_.find(dirPath);
  if (it != subdirs_.end()) {
    return it->second;
  }
  // Build nested directories one segment at a time, caching each level so two
  // histograms in the same subdirectory do not create it twice.
  const std::size_t slash = dirPath.rfind('/');
  art::TFileDirectory& parent =
      (slash == std::string::npos) ? *dir_ : directoryFor(dirPath.substr(0, slash));
  const std::string leaf =
      (slash == std::string::npos) ? dirPath : dirPath.substr(slash + 1);
  auto inserted = subdirs_.emplace(dirPath, parent.mkdir(leaf));
  return inserted.first->second;
}

const DQMSegmentation::Rule& DQMSegmentation::ruleFor(const std::string& path) const
{
  for (const auto& rule : config_.rules) {
    if (globMatch(rule.match, path)) {
      return rule;
    }
  }
  return defaultRule_;
}

TH1* DQMSegmentation::create(Entry& entry, const std::string& dirPath,
                             const std::string& name, bool persist)
{
  if (persist) {
    return entry.makeInDir(directoryFor(dirPath), name, entry.title);
  }
  // Belongs to no TDirectory, so nothing writes or deletes it behind our back.
  auto h = entry.makeOwned(name, entry.title);
  TH1* raw = h.get();
  owned_.push_back(std::move(h));
  return raw;
}

void DQMSegmentation::dropOwned(TH1* h)
{
  auto it = std::find_if(owned_.begin(), owned_.end(),
                         [h](const std::unique_ptr<TH1>& p) { return p.get() == h; });
  if (it != owned_.end()) {
    owned_.erase(it);
  }
}

void DQMSegmentation::refreshTargets(Entry& entry)
{
  entry.targets.clear();
  entry.targets.add(entry.job);
  entry.targets.add(entry.subrunLive);
  entry.targets.add(entry.windowLive);
  if (!entry.subBlocks.empty()) {
    entry.targets.add(entry.subBlocks[entry.cur]);
  }
}

DQMSegmentation::Entry* DQMSegmentation::addEntry(const std::string& path,
                                                  const std::string& title,
                                                  DirFactory makeInDir,
                                                  OwnedFactory makeOwned)
{
  if (!dir_) {
    return nullptr;
  }
  auto existing = byPath_.find(path);
  if (existing != byPath_.end()) {
    return existing->second;
  }

  const Rule& rule = ruleFor(path);
  if (!rule.enabled) {
    return nullptr;  //handle stays empty; every fill through it is a no-op
  }

  auto owner = std::make_unique<Entry>();
  Entry& entry = *owner;
  entry.path = path;
  entry.title = title;
  entry.rule = rule;
  entry.makeInDir = std::move(makeInDir);
  entry.makeOwned = std::move(makeOwned);

  const std::size_t slash = path.rfind('/');
  entry.dirPath = (slash == std::string::npos) ? "" : path.substr(0, slash);
  entry.base = (slash == std::string::npos) ? path : path.substr(slash + 1);

  const std::string segmentPath =
      entry.dirPath.empty() ? config_.segmentDir
                            : entry.dirPath + "/" + config_.segmentDir;

  if (rule.job) {
    entry.job = create(entry, entry.dirPath, entry.base, rule.jobPersist);
  }

  if (rule.subrun.enabled) {
    entry.subrunLive = create(entry, entry.dirPath, entry.base + "_sub",
                              rule.subrun.persistLive);
  }

  if (rule.window.enabled) {
    const std::string liveName =
        rule.liveName.empty() ? entry.base + "_last" : rule.liveName;
    entry.windowLive =
        create(entry, entry.dirPath, liveName, rule.window.persistLive);

    const int nSub = std::max(rule.window.subdivisions, 1);
    entry.subBlocks.reserve(nSub);
    entry.subBlockRange.resize(nSub);
    for (int i = 0; i < nSub; ++i) {
      // Internal accumulators: never written, so their names only have to be
      // unique within the process.
      entry.subBlocks.push_back(
          create(entry, entry.dirPath, entry.base + Form("_subblock%d", i), false));
    }
    entry.subBlockWidth =
        std::max<unsigned long long>(rule.window.span / static_cast<unsigned long long>(nSub), 1);

    const int nKeep = std::max(rule.window.keep, 0);
    entry.windowArchive.reserve(nKeep);
    entry.archiveRange.resize(nKeep);
    for (int i = 0; i < nKeep; ++i) {
      entry.windowArchive.push_back(create(entry, segmentPath,
                                           entry.base + Form("_prev%d", i + 1),
                                           rule.window.persist));
    }
  }

  refreshTargets(entry);
  labelAll(entry);

  Entry* raw = owner.get();
  entries_.push_back(std::move(owner));
  byPath_[path] = raw;
  if (raw->subrunLive != nullptr || raw->windowLive != nullptr) {
    segmented_.push_back(raw);
  }
  return raw;
}

unsigned long long DQMSegmentation::clockFor(const Entry& entry) const
{
  switch (entry.rule.window.unit) {
    case Unit::Event:
      return static_cast<unsigned long long>(currentEvent_);
    case Unit::Ewt:
      return haveEwt_ ? currentEwt_ : static_cast<unsigned long long>(currentEvent_);
    case Unit::SubRun:
      return subrunClock_;
  }
  return static_cast<unsigned long long>(currentEvent_);
}

void DQMSegmentation::Advance(std::size_t eventIndex, std::optional<uint64_t> ewt)
{
  currentEvent_ = eventIndex;
  if (ewt.has_value()) {
    currentEwt_ = *ewt;
    haveEwt_ = true;
  }

  for (Entry* owner : segmented_) {
    Entry& entry = *owner;
    const unsigned long long clock = clockFor(entry);

    if (entry.rule.window.enabled && entry.rule.window.unit == Unit::Ewt &&
        !haveEwt_ && !ewtFallback_) {
      ewtFallback_ = true;
      mf::LogWarning("DQMSegmentation")
          << "a window rule asks for unit \"ewt\" but the input carries no event "
          << "window tag (typical MC). Falling back to counting events for the "
          << "window clock. Reported once per job.";
    }

    if (!entry.subBlocks.empty()) {
      advanceWindow(entry, clock);
    }

    note(entry.jobRange, eventIndex, clock, run_, subrun_);
    if (entry.subrunLive != nullptr) {
      note(entry.subrunRange, eventIndex, clock, run_, subrun_);
    }
    if (entry.windowLive != nullptr) {
      note(entry.spanRange, eventIndex, clock, run_, subrun_);
      note(entry.subBlockRange[entry.cur], eventIndex, clock, run_, subrun_);
      // Extended here rather than only at rotations, so the live copy's title
      // is never behind the data in it. rebuildLive recomputes it from the
      // ring when a sub-block is evicted.
      note(entry.liveRange, eventIndex, clock, run_, subrun_);
    }
  }
}

void DQMSegmentation::advanceWindow(Entry& entry, unsigned long long clock)
{
  if (!entry.windowStarted) {
    entry.windowStarted = true;
    entry.subBlockStartClock = clock;
    return;
  }
  if (clock < entry.subBlockStartClock) {
    // The clock went backwards (a new run, or an EWT reset). Everything in the
    // ring belongs to a different regime; start over rather than mislabel it.
    for (std::size_t i = 0; i < entry.subBlocks.size(); ++i) {
      entry.subBlocks[i]->Reset("ICES");
      entry.subBlockRange[i] = Range{};
    }
    entry.spanRange = Range{};
    entry.cur = 0;
    entry.completedSubBlocks = 0;
    entry.subBlockStartClock = clock;
    rebuildLive(entry);
    refreshTargets(entry);
    return;
  }

  const unsigned long long elapsed = clock - entry.subBlockStartClock;
  std::size_t nRotations = static_cast<std::size_t>(elapsed / entry.subBlockWidth);
  if (nRotations == 0) {
    return;
  }
  // More than a full ring of silence means nothing retained is still in the
  // window; rotating the whole ring once is both correct and bounded.
  const std::size_t maxRotations = entry.subBlocks.size();
  const bool capped = nRotations > maxRotations;
  if (capped) {
    nRotations = maxRotations;
  }
  for (std::size_t i = 0; i < nRotations; ++i) {
    rotateSubBlock(entry, clock);
  }
  entry.subBlockStartClock =
      capped ? clock
             : entry.subBlockStartClock +
                   nRotations * entry.subBlockWidth;
}

void DQMSegmentation::rotateSubBlock(Entry& entry, unsigned long long /*clock*/)
{
  ++entry.completedSubBlocks;
  if (entry.completedSubBlocks % entry.subBlocks.size() == 0) {
    // The ring now holds exactly one span of completed sub-blocks, and the live
    // copy is their sum, so it is the span to archive.
    archiveSpan(entry);
    entry.spanRange = Range{};
  }

  entry.cur = (entry.cur + 1) % entry.subBlocks.size();
  entry.subBlocks[entry.cur]->Reset("ICES");
  entry.subBlockRange[entry.cur] = Range{};
  rebuildLive(entry);
  refreshTargets(entry);
}

void DQMSegmentation::archiveSpan(Entry& entry)
{
  if (entry.windowArchive.empty()) {
    return;
  }
  Range finished = entry.spanRange;
  finished.complete = true;

  // Contents shift down the fixed _prevN names so the online GUI can subscribe
  // to a stable name; _prev1 is always the most recent completed span.
  for (std::size_t j = entry.windowArchive.size() - 1; j > 0; --j) {
    entry.windowArchive[j]->Reset("ICES");
    entry.windowArchive[j]->Add(entry.windowArchive[j - 1]);
    entry.archiveRange[j] = entry.archiveRange[j - 1];
    label(entry, entry.windowArchive[j], entry.archiveRange[j], "window",
          static_cast<int>(j + 1));
  }
  entry.windowArchive[0]->Reset("ICES");
  if (entry.windowLive != nullptr) {
    entry.windowArchive[0]->Add(entry.windowLive);
  }
  entry.archiveRange[0] = finished;
  label(entry, entry.windowArchive[0], entry.archiveRange[0], "window", 1);
}

void DQMSegmentation::rebuildLive(Entry& entry)
{
  if (entry.windowLive == nullptr) {
    return;
  }
  // Re-summing rather than subtracting the evicted block keeps the bin errors
  // and the entry count correct.
  entry.windowLive->Reset("ICES");
  entry.liveRange = Range{};
  for (std::size_t i = 0; i < entry.subBlocks.size(); ++i) {
    entry.windowLive->Add(entry.subBlocks[i]);
    const Range& r = entry.subBlockRange[i];
    if (!r.started) {
      continue;
    }
    if (!entry.liveRange.started) {
      entry.liveRange = r;
      continue;
    }
    entry.liveRange.firstEvent = std::min(entry.liveRange.firstEvent, r.firstEvent);
    entry.liveRange.lastEvent = std::max(entry.liveRange.lastEvent, r.lastEvent);
    entry.liveRange.firstClock = std::min(entry.liveRange.firstClock, r.firstClock);
    entry.liveRange.lastClock = std::max(entry.liveRange.lastClock, r.lastClock);
    entry.liveRange.nEvents += r.nEvents;
    entry.liveRange.subrun = r.subrun;
    entry.liveRange.run = r.run;
  }
  label(entry, entry.windowLive, entry.liveRange, "windowLive", 0);
}

void DQMSegmentation::BeginSubRun(int run, int subrun)
{
  run_ = run;
  subrun_ = subrun;
}

void DQMSegmentation::EndSubRun()
{
  for (Entry* owner : segmented_) {
    Entry& entry = *owner;
    if (entry.subrunLive == nullptr) {
      continue;
    }
    Range finished = entry.subrunRange;
    finished.complete = true;
    // Name and label from the registry's current run/subrun rather than what
    // the first event of the range carried, so a caller that only calls
    // EndSubRun still gets an unambiguous name.
    finished.run = run_;
    finished.subrun = subrun_;

    if (finished.started) {
      const std::string archiveDir =
          entry.dirPath.empty() ? config_.subrunDir
                                : entry.dirPath + "/" + config_.subrunDir;
      const std::string name =
          entry.base + Form("_r%06d_s%06d", finished.run, finished.subrun);
      TH1* archive = create(entry, archiveDir, name, entry.rule.subrun.persist);
      archive->Reset("ICES");
      archive->Add(entry.subrunLive);
      label(entry, archive, finished, "subrun",
            static_cast<int>(entry.nSubRunsArchived + 1));
      ++entry.nSubRunsArchived;

      // `keep` bounds the in-memory history. A copy already written to the file
      // is what persisting asked for, so it is never dropped.
      if (!entry.rule.subrun.persist && entry.rule.subrun.keep >= 0) {
        entry.subrunArchive.push_back(archive);
        while (entry.subrunArchive.size() >
               static_cast<std::size_t>(entry.rule.subrun.keep)) {
          dropOwned(entry.subrunArchive.front());
          entry.subrunArchive.erase(entry.subrunArchive.begin());
        }
      } else {
        entry.subrunArchive.push_back(archive);
      }
    }

    entry.subrunLive->Reset("ICES");
    entry.subrunRange = Range{};
    label(entry, entry.subrunLive, entry.subrunRange, "subrunLive", 0);
  }
  ++subrunClock_;
}

void DQMSegmentation::label(Entry& entry, TH1* h, const Range& range,
                            const char* mode, int index)
{
  if (h == nullptr) {
    return;
  }
  // A histogram with only a job copy has no sibling to be confused with, and
  // annotating it would change the output of an unconfigured job. Leave it
  // exactly as the helper booked it.
  if (entry.subrunLive == nullptr && entry.windowLive == nullptr) {
    return;
  }

  std::string tag;
  if (config_.annotateTitles) {
    const char* unit = unitName(entry.rule.window.unit);
    std::ostringstream os;
    if (std::string(mode) == "job") {
      os << "[job: events " << range.firstEvent << "-" << range.lastEvent << ", "
         << range.nEvents << " events]";
    } else if (std::string(mode) == "subrunLive") {
      os << "[run " << range.run << " subrun " << range.subrun
         << " in progress: " << range.nEvents << " events]";
    } else if (std::string(mode) == "subrun") {
      os << "[run " << range.run << " subrun " << range.subrun << ": "
         << range.nEvents << " events]";
    } else if (std::string(mode) == "windowLive") {
      os << "[last " << entry.rule.window.span << " " << unit << ": "
         << range.firstClock << "-" << range.lastClock << ", " << range.nEvents
         << " events]";
    } else {
      os << "[previous span " << index << " of " << entry.windowArchive.size()
         << ": " << unit << " " << range.firstClock << "-" << range.lastClock
         << ", " << range.nEvents << " events]";
    }
    if (!range.started) {
      os.str("");
      os << "[" << mode << ": empty]";
    }
    tag = os.str();
  }
  h->SetTitle(annotate(entry.title, tag).c_str());

  if (!config_.stampMetadata) {
    return;
  }
  TList* functions = h->GetListOfFunctions();
  if (functions == nullptr) {
    return;
  }
  // Reset("ICES") clears the function list along with the contents and stats,
  // so every reset is followed by a call here that puts the stamp back.
  if (TObject* old = functions->FindObject("dqmSegment")) {
    functions->Remove(old);
    delete old;
  }
  functions->Add(new TNamed("dqmSegment", stampText(range, mode, index).c_str()));
}

void DQMSegmentation::labelAll(Entry& entry)
{
  label(entry, entry.job, entry.jobRange, "job", 0);
  label(entry, entry.subrunLive, entry.subrunRange, "subrunLive", 0);
  label(entry, entry.windowLive, entry.liveRange, "windowLive", 0);
  for (std::size_t j = 0; j < entry.windowArchive.size(); ++j) {
    label(entry, entry.windowArchive[j], entry.archiveRange[j], "window",
          static_cast<int>(j + 1));
  }
}

void DQMSegmentation::RefreshLabels()
{
  for (Entry* owner : segmented_) {
    Entry& entry = *owner;
    label(entry, entry.job, entry.jobRange, "job", 0);
    label(entry, entry.subrunLive, entry.subrunRange, "subrunLive", 0);
    label(entry, entry.windowLive, entry.liveRange, "windowLive", 0);
  }
}

void DQMSegmentation::Finalize()
{
  for (Entry* owner : segmented_) {
    Entry& entry = *owner;
    entry.jobRange.complete = true;
    entry.liveRange.complete = true;
    labelAll(entry);
  }
}

std::vector<TH1*> DQMSegmentation::copies(const std::string& path) const
{
  std::vector<TH1*> out;
  auto it = byPath_.find(path);
  if (it == byPath_.end()) {
    return out;
  }
  const Entry& entry = *it->second;
  if (entry.job != nullptr) {
    out.push_back(entry.job);
  }
  if (entry.subrunLive != nullptr) {
    out.push_back(entry.subrunLive);
  }
  for (TH1* h : entry.subrunArchive) {
    out.push_back(h);
  }
  if (entry.windowLive != nullptr) {
    out.push_back(entry.windowLive);
  }
  for (TH1* h : entry.windowArchive) {
    out.push_back(h);
  }
  return out;
}

std::vector<TH1*> DQMSegmentation::allCopies() const
{
  std::vector<TH1*> out;
  for (const auto& owner : entries_) {
    for (TH1* h : copies(owner->path)) {
      out.push_back(h);
    }
  }
  return out;
}

std::map<std::string, std::vector<TH1*>> DQMSegmentation::publishedCopies() const
{
  std::map<std::string, std::vector<TH1*>> out;
  for (const auto& owner : entries_) {
    if (!owner->rule.publish) {
      continue;
    }
    for (TH1* h : copies(owner->path)) {
      // No group: the copy is published under its own name, so the online
      // names stay the ones a GUI already subscribes to.
      out[owner->rule.group.empty() ? h->GetName() : owner->rule.group].push_back(h);
    }
  }
  return out;
}

TH1* DQMSegmentation::live(const std::string& path) const
{
  auto it = byPath_.find(path);
  if (it == byPath_.end()) {
    return nullptr;
  }
  return it->second->windowLive;
}

} // namespace mu2e
