// Ed Callaghan
// Container for CaloDigiWrappers, with functionality to resolve overlapping digitization windows
// September 2025

#include "Offline/CaloMC/inc/CaloDigiWrapperCollection.hh"

namespace mu2e{
  // forward size query to underlying container
  size_t CaloDigiWrapperCollection::size() const{
    auto rv = _wrappers.size();
    return rv;
  }

  // forward iterator access to underlying wrappers
  CDWC_iterator CaloDigiWrapperCollection::begin(){
    auto rv = _wrappers.begin();
    return rv;
  }

  CDWC_const_iterator CaloDigiWrapperCollection::begin() const{
    auto rv = _wrappers.begin();
    return rv;
  }

  CDWC_iterator CaloDigiWrapperCollection::end(){
    auto rv = _wrappers.end();
    return rv;
  }

  CDWC_const_iterator CaloDigiWrapperCollection::end() const{
    auto rv = _wrappers.end();
    return rv;
  }

  // forward lookups to underlying container
  CaloDigiWrapper& CaloDigiWrapperCollection::operator[](size_t i){
    auto& rv = _wrappers[i];
    return rv;
  }

  const CaloDigiWrapper& CaloDigiWrapperCollection::operator[](size_t i) const{
    const auto& rv = _wrappers[i];
    return rv;
  }

  CaloDigiWrapper& CaloDigiWrapperCollection::front(){
    auto& rv = _wrappers.front();
    return rv;
  }

  const CaloDigiWrapper& CaloDigiWrapperCollection::front() const{
    const auto& rv = _wrappers.front();
    return rv;
  }

  CaloDigiWrapper& CaloDigiWrapperCollection::back(){
    auto& rv = _wrappers.back();
    return rv;
  }

  const CaloDigiWrapper& CaloDigiWrapperCollection::back() const{
    const auto& rv = _wrappers.back();
    return rv;
  }

  // insertion
  void CaloDigiWrapperCollection::Append(const CaloDigiCollection& digis){
    for (size_t i = 0 ; i < digis.size() ; i++){
      const auto& digi = digis.at(i);
      CaloDigiWrapper wrapper(digi);
      _wrappers.emplace_back(digi);
      continue;
    }
  }

  void CaloDigiWrapperCollection::Append(const CaloDigiWrapper& wrapper){
    _wrappers.push_back(wrapper);
  }

  // accessor
  std::unique_ptr<CaloDigiCollection> CaloDigiWrapperCollection::GetDigis() const{
    auto rv = std::make_unique<CaloDigiCollection>();
    rv->reserve(_wrappers.size());
    for (const auto& wrapper: _wrappers){
      const auto& digi = wrapper.Digi();
      rv->push_back(digi);
    }

    return rv;
  }

  // collision resolution
  // auxiliary comparison to facilitate time-sorting
  bool compare_t0s(const CaloDigiWrapper* lhs, const CaloDigiWrapper* rhs){
    bool rv = (lhs->Digi().t0() < rhs->Digi().t0());
    return rv;
  }

  // ejc: this overlaps _greatly_ with CaloDigiWrapperCollection
  // TODO this assumes a fixed overlap window, which does not apply for calo
  // here we need variable-length chain links --- -.-
  void CaloDigiWrapperCollection::ResolveCollisions(CaloDigiWrapperCollection& rv){
    // identify time-overlapped chains: this is a 3 step process
    // first, partition wrappers according to SiPMID_t
    std::map<SiPMID_t, CaloDigiWrapperCollection> unsorted_map;
    for (const auto& wrapper: (*this)){
      const auto& digi = wrapper.Digi();
      SiPMID_t id = digi.SiPMID();
      unsorted_map[id].Append(wrapper);
    }

    // next, (approximately) sort all subsets via time, to avoid
    // pathological situations where timing buckets are misconstructed
    // we shamefully and begrudgingly do this indirectly via bare pointers
    std::map<SiPMID_t, CaloDigiWrapperCollection> wrappers_map;
    for (const auto& pair: unsorted_map){
      const auto& id = pair.first;
      auto& wrappers = pair.second;
      std::vector<CaloDigiWrapper*> sortable(wrappers.size());
      for (size_t i = 0 ; i < wrappers.size() ; i++){
        sortable[i] = const_cast<CaloDigiWrapper*>(&wrappers[i]);
      }
      std::sort(sortable.begin(), sortable.end(), compare_t0s);
      CaloDigiWrapperCollection sorted;
      for (size_t i = 0 ; i < wrappers.size() ; i++){
        sorted.Append(*sortable[i]);
      }
      wrappers_map.emplace(id, sorted);
    }

    // then, filter each subset into buckets based on overlapping time windows
    std::map< SiPMID_t, TimeBasedBuckets<CaloDigiWrapper> > buckets_map;
    for (const auto& pair: wrappers_map){
      const auto& id = pair.first;
      const auto& wrappers = pair.second;
      TimeBasedBuckets<CaloDigiWrapper> buckets;
      for (const auto& wrapper: wrappers){
        double window = static_cast<double>(wrapper.Digi().waveform().size());
        buckets.Insert(wrapper, window);
      }
      buckets_map.emplace(id, buckets);
    }

    // finally, recast each bucket as a regular container
    // and resolve the set of degenerate digis into one
    for (const auto& pair: buckets_map){
      const auto& buckets = pair.second;
      for (const auto& bucket: buckets){
        CaloDigiWrapperCollection wrappers;
        for (const auto& wrapper: bucket){
          wrappers.Append(wrapper);
        }
        this->ResolveCollision(wrappers, rv);
      }
    }
  }

  void CaloDigiWrapperCollection::ResolveCollision(CaloDigiWrapperCollection& collided, CaloDigiWrapperCollection& rv){
    // if only one digi present, then nothing to do
    if (collided.size() < 2){
      rv.Append(collided.front());
      return;
    }

    // multiple digis must be summed into a single waveform
    // here we determine the length of that waveform, by finding
    // how far past the end of the first it extends, if at all
    const auto& first = collided.front().Digi();
    size_t length = 0;
    for (const auto& wrapper: collided){
      const auto& digi = wrapper.Digi();
      size_t proposed = (digi.t0() - first.t0()) + digi.waveform().size();
      if (length < proposed){
        length = proposed;
      }
    }
    std::vector<sample_t> samples(length, 0);

    // sum each individual waveform into the total
    // TODO this does not account for saturation
    for (const auto& wrapper: collided){
      const auto& digi = wrapper.Digi();
      const auto& waveform = digi.waveform();
      const size_t shift = digi.t0() - first.t0();
      for (size_t i = 0 ; i < waveform.size() ; i++){
        samples[i + shift] += waveform[i];
      }
    }

    // now find the peak position
    pos_t peakpos = 0;
    sample_t peak = 0;
    for (size_t i = 0 ; i < samples.size() ; i++){
      const auto sample = samples[i];
      if (peak < sample){
        peakpos = i;
        peak = sample;
      }
    }

    const auto id = first.SiPMID();
    const auto t0 = first.t0();
    const auto digi = CaloDigi(id, t0, samples, peakpos);
    const auto wrapper = CaloDigiWrapper(digi);
    rv.Append(wrapper);
  }
}
