// FamilySplitSummary: per-SubRun counts of a cosmic stage-2 job made with Geant4 secondary splitting, taken by
// FamilySplitCounter before any event selection that follows it. FamilySplitSource (with countsLabel) weights
// the stripped events with them instead of with the events it finds in its input file, so that stage 2 may
// write only the events a later selection needs:
//   stripped-event weight = 1 - sumWeights / nKept
// nEvents: events the counter saw; nKept: those kept by the source's muon-decay rule; sumWeights: the sum of
// weights (1/N) of the forced families of the kept events.
#ifndef MCDataProducts_inc_FamilySplitSummary_hh
#define MCDataProducts_inc_FamilySplitSummary_hh

namespace mu2e {

  class FamilySplitSummary {
  public:
    FamilySplitSummary() = default;
    FamilySplitSummary(unsigned long nEvents, unsigned long nKept, double sumWeights)
      : nEvents_(nEvents), nKept_(nKept), sumWeights_(sumWeights) {}

    unsigned long nEvents() const { return nEvents_; }
    unsigned long nKept() const { return nKept_; }
    double sumWeights() const { return sumWeights_; }

    // art calls this when it combines the same SubRun from several fragments
    void aggregate(FamilySplitSummary const& other) {
      nEvents_ += other.nEvents_;
      nKept_ += other.nKept_;
      sumWeights_ += other.sumWeights_;
    }

  private:
    unsigned long nEvents_ = 0;
    unsigned long nKept_ = 0;
    double sumWeights_ = 0.;
  };

} // namespace mu2e

#endif /* MCDataProducts_inc_FamilySplitSummary_hh */
