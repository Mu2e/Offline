#ifndef DQMHelpers_inc_DQMHist_hh
#define DQMHelpers_inc_DQMHist_hh
// Lightweight handles returned by DQMSegmentation::book*. A handle fans one
// Fill out to every copy of that histogram the FHiCL segmentation asked for
// (job, current subrun, rolling window, current window sub-block), so helper
// code fills a histogram the same way whether it is segmented or not.
//
// Original Author: R. Mina

#include "TH1.h"

#include <functional>

namespace mu2e {

// Fill targets for one booked histogram. DQMSegmentation owns this and
// rewrites h[] when the window ring rotates; handles read it at fill time, so
// they never go stale. Four is the most copies a single fill can reach.
struct DQMHistTargets {
  TH1* h[4]{nullptr, nullptr, nullptr, nullptr};
  int n{0};

  void clear() { h[0] = h[1] = h[2] = h[3] = nullptr; n = 0; }
  void add(TH1* p)
  {
    if (p != nullptr && n < 4) {
      h[n++] = p;
    }
  }
};

// Common part of both handles. `H*` conversion and operator-> both give the
// job copy, so existing accessors and null checks keep working unchanged.
template <class H>
class DQMHistBase {
public:
  DQMHistBase() = default;
  explicit DQMHistBase(const DQMHistTargets* targets, H* job) :
      targets_(targets), job_(job)
  {}

  H* get() const { return job_; }
  operator H*() const { return job_; }
  H* operator->() const { return job_; }

  // Apply styling at book time to every copy, not just the one that is written.
  void ForEach(const std::function<void(H*)>& f) const
  {
    if (targets_ == nullptr) {
      return;
    }
    for (int i = 0; i < targets_->n; ++i) {
      f(static_cast<H*>(targets_->h[i]));
    }
  }

  void AddBinContent(int bin, double w) const
  {
    if (targets_ == nullptr) {
      return;
    }
    for (int i = 0; i < targets_->n; ++i) {
      targets_->h[i]->AddBinContent(bin, w);
    }
  }

protected:
  const DQMHistTargets* targets_{nullptr};
  H* job_{nullptr};
};

// 1D: Fill(x) and Fill(x, weight).
template <class H>
class DQMHist1 : public DQMHistBase<H> {
public:
  using DQMHistBase<H>::DQMHistBase;

  void Fill(double x) const
  {
    const DQMHistTargets* t = this->targets_;
    if (t == nullptr) {
      return;
    }
    for (int i = 0; i < t->n; ++i) {
      t->h[i]->Fill(x);
    }
  }

  void Fill(double x, double w) const
  {
    const DQMHistTargets* t = this->targets_;
    if (t == nullptr) {
      return;
    }
    for (int i = 0; i < t->n; ++i) {
      t->h[i]->Fill(x, w);
    }
  }
};

// 2D: Fill(x, y) and Fill(x, y, weight). Kept a separate type so a 2-argument
// fill cannot silently mean "x with weight" on one hist and "x, y" on another.
template <class H>
class DQMHist2 : public DQMHistBase<H> {
public:
  using DQMHistBase<H>::DQMHistBase;

  void Fill(double x, double y) const
  {
    const DQMHistTargets* t = this->targets_;
    if (t == nullptr) {
      return;
    }
    for (int i = 0; i < t->n; ++i) {
      static_cast<H*>(t->h[i])->Fill(x, y);
    }
  }

  void Fill(double x, double y, double w) const
  {
    const DQMHistTargets* t = this->targets_;
    if (t == nullptr) {
      return;
    }
    for (int i = 0; i < t->n; ++i) {
      static_cast<H*>(t->h[i])->Fill(x, y, w);
    }
  }
};

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMHist_hh */
