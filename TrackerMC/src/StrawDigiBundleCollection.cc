// Containerize StrawDigiBundles, with a extra functionality for resolving
// situations where multiple digis have overlapping digitization windows
//
// Ed Callaghan, 2023

#include <Offline/TrackerMC/inc/StrawDigiBundleCollection.hh>
#include <Offline/TrackerMC/inc/TimeBasedBucket.hh>

namespace mu2e{
  // convenience methods for accepting new StrawDigiBundles
  // from different source situations
  void StrawDigiBundleCollection::Append(const StrawDigiBundle bundle){
    this->push_back(bundle);
  }

  // utility functions
  void StrawDigiBundleCollection::FillStrawDigis(StrawDigiCollection& rv){
    rv.resize(this->size());
    for (size_t i = 0 ; i < this->size() ; i++){
      rv[i] = this->at(i).GetStrawDigi();
    }
  }

  void StrawDigiBundleCollection::FillStrawDigiADCWaveforms(StrawDigiADCWaveformCollection& rv){
    rv.resize(this->size());
    for (size_t i = 0 ; i < this->size() ; i++){
      rv[i] = this->at(i).GetStrawDigiADCWaveform();
    }
  }

  void StrawDigiBundleCollection::FillStrawDigiMCs(StrawDigiMCCollection& rv){
    rv.resize(this->size());
    for (size_t i = 0 ; i < this->size() ; i++){
      rv[i] = this->at(i).GetStrawDigiMC();
    }
  }

  StrawDigiBundleCollection StrawDigiBundleCollection::operator+= (const StrawDigiBundleCollection& other){
    for (auto bundle: other){
      this->push_back(bundle);
    }
    return (*this);
  }

  // accessors, as objects
  StrawDigiCollection StrawDigiBundleCollection::GetStrawDigis(){
    StrawDigiCollection rv;
    this->FillStrawDigis(rv);
    return rv;
  }

  StrawDigiADCWaveformCollection StrawDigiBundleCollection::GetStrawDigiADCWaveforms(){
    StrawDigiADCWaveformCollection rv;
    this->FillStrawDigiADCWaveforms(rv);
    return rv;
  }

  StrawDigiMCCollection StrawDigiBundleCollection::GetStrawDigiMCs(){
    StrawDigiMCCollection rv;
    this->FillStrawDigiMCs(rv);
    return rv;
  }

  // accessors, as pointers
  std::unique_ptr<StrawDigiCollection> StrawDigiBundleCollection::GetStrawDigiPtrs(){
    auto rv = std::make_unique<StrawDigiCollection>();
    this->FillStrawDigis(*rv);
    return rv;
  }

  std::unique_ptr<StrawDigiADCWaveformCollection> StrawDigiBundleCollection::GetStrawDigiADCWaveformPtrs(){
    auto rv = std::make_unique<StrawDigiADCWaveformCollection>();
    this->FillStrawDigiADCWaveforms(*rv);
    return rv;
  }

  std::unique_ptr<StrawDigiMCCollection> StrawDigiBundleCollection::GetStrawDigiMCPtrs(){
    auto rv = std::make_unique<StrawDigiMCCollection>();
    this->FillStrawDigiMCs(*rv);
    return rv;
  }

  // auxiliary comparison to facilitate time-sorting
  bool compare_tdcs(const StrawDigiBundle* lhs, const StrawDigiBundle* rhs){
    bool rv = (lhs->GetStrawDigi().TDC() < rhs->GetStrawDigi().TDC());
    return rv;
  }

  StrawDigiBundleCollection StrawDigiBundleCollection::ResolveCollisions(const StrawElectronics& conditions){
    // identify time-overlapped chains: this is a 3 step process
    // first, partition bundles according to StrawId
    std::map<StrawId, StrawDigiBundleCollection> unsorted_map;
    for (const auto& bundle: (*this)){
      const auto& digi = bundle.GetStrawDigi();
      StrawId id = digi.strawId();
      unsorted_map[id].push_back(bundle);
    }

    // next, (approximately) sort all subsets via time, to avoid
    // pathological situations where timing buckets are misconstructed
    // we shamefully and begrudgingly do this indirectly via bare pointers
    std::map<StrawId, StrawDigiBundleCollection> bundles_map;
    for (const auto& pair: unsorted_map){
      const auto& id = pair.first;
      auto& bundles = pair.second;
      std::vector<StrawDigiBundle*> sortable(bundles.size());
      for (size_t i = 0 ; i < bundles.size() ; i++){
        sortable[i] = const_cast<StrawDigiBundle*>(&bundles[i]);
      }
      std::sort(sortable.begin(), sortable.end(), compare_tdcs);
      std::vector<StrawDigiBundle> sorted;
      for (size_t i = 0 ; i < bundles.size() ; i++){
        sorted.push_back(*sortable[i]);
      }
      bundles_map.emplace(id, sorted);
    }

    // then, filter each subset into buckets based on overlapping time windows
    double window = conditions.deadTimeDigital() / conditions.tdcLSB();
    std::map< StrawId, TimeBasedBuckets<StrawDigiBundle> > buckets_map;
    for (const auto& pair: bundles_map){
      const auto& id = pair.first;
      const auto& bundles = pair.second;
      TimeBasedBuckets<StrawDigiBundle> buckets(window);
      for (const auto& bundle: bundles){
        buckets.Insert(bundle);
      }
      buckets_map.emplace(id, buckets);
    }

    // finally, recast each bucket as a regular container
    // and resolve the set of degenerate digis into one
    StrawDigiBundleCollection rv;
    for (const auto& pair: buckets_map){
      const auto& buckets = pair.second;
      for (const auto& bucket: buckets){
        StrawDigiBundleCollection bundles;
        for (const auto& bundle: bucket){
          bundles.Append(bundle);
        }
        auto resolved = this->ResolveCollision(bundles);
        rv += resolved;
      }
    }

    return rv;
  }

  StrawDigiBundleCollection StrawDigiBundleCollection::ResolveCollision(StrawDigiBundleCollection& collided){
    StrawDigiBundleCollection rv;
    auto first = collided[0];
    // update StrawDigiMC component to reflect that the MC information may be tainted or incomplete
    auto mc = StrawDigiMC(first.GetStrawDigiMC(), StrawDigiMC::PartiallyValid);
    auto updated = StrawDigiBundle(first.GetStrawDigi(), first.GetStrawDigiADCWaveform(), mc);
    rv.push_back(updated);
    return rv;
  }
}
