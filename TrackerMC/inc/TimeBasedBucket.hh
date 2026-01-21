// Ed Callaghan
// Time-based bucketing, of e.g. StrawGasSteps and StrawDigis into
// groups of overlapping quasi-digitization windows
// February 2024

// stl
#include <vector>

namespace mu2e{
  // alias iterators of underlying container, for easy loops
  template<typename T>
  using TBB_iterator = std::vector<T>::iterator;
  template<typename T>
  using TBB_const_iterator = std::vector<T>::const_iterator;

  // container-type which is aware of whether a candidate item should be
  // added, or not
  template<typename T>
  class TimeBasedBucket{
    public:
      TimeBasedBucket();

      // forwarded calls to underlying container
      size_t size() const;
      TBB_iterator<T> begin();
      TBB_const_iterator<T> begin() const;
      TBB_iterator<T> end();
      TBB_const_iterator<T> end() const;

      // whether a candidate item "fits" into this bucket's gross time window
      bool Accepts(T, double);
      // insertion --- assumes that appends are pre-sorted
      void Append(T, double);
    protected:
      std::vector<T> _container;
      std::vector<double> _windows;
    private:
      /**/
  };

  // alias iterators of underlying container, for easy loops
  template<typename T>
  using TBBS_iterator = std::vector< TimeBasedBucket<T> >::iterator;
  template<typename T>
  using TBBS_const_iterator = std::vector< TimeBasedBucket<T> >::const_iterator;

  // a container of such containers, which ensures that every item considered
  // either lands in an existing bucket, or defines its own
  template <typename T>
  class TimeBasedBuckets{
    public:
      TimeBasedBuckets();

      // forwarded calls to underlying container
      size_t size() const;
      TBBS_iterator<T> begin();
      TBBS_const_iterator<T> begin() const;
      TBBS_iterator<T> end();
      TBBS_const_iterator<T> end() const;

      // insertion
      void Insert(T, double);
    protected:
      std::vector< TimeBasedBucket<T> > _buckets;
    private:
      /**/
  };

  template <typename T>
  TimeBasedBucket<T>::TimeBasedBucket(){
    /**/
  }

  // forward size query to underlying container
  template <typename T>
  size_t TimeBasedBucket<T>::size() const{
    auto rv = _container.size();
    return rv;
  }

  // forward iterator access to underlying container
  template <typename T>
  TBB_iterator<T> TimeBasedBucket<T>::begin(){
    auto rv = _container.begin();
    return rv;
  }

  template <typename T>
  TBB_const_iterator<T> TimeBasedBucket<T>::begin() const{
    auto rv = _container.begin();
    return rv;
  }

  template <typename T>
  TBB_iterator<T> TimeBasedBucket<T>::end(){
    auto rv = _container.end();
    return rv;
  }

  template <typename T>
  TBB_const_iterator<T> TimeBasedBucket<T>::end() const{
    auto rv = _container.end();
    return rv;
  }

  // whether a candidate item "fits" into this bucket's gross time window
  template <typename T>
  bool TimeBasedBucket<T>::Accepts(T candidate, double window){
    // empty bucket rejects nothing
    if (this->size() < 1){
      return true;
    }

    // otherwise, check between first/last times, accounting for windowing
    double lower = _container.front().time() - window;
    double upper = _container.back().time() + _windows.back();
    double now = candidate.time();
    bool rv = (lower <= now) && (now <= upper);
    return rv;
  }

  // insertion --- assumes that appends are pre-sorted
  template<typename T>
  void TimeBasedBucket<T>::Append(T candidate, double window){
    _container.push_back(candidate);
    _windows.push_back(window);
  }

  template<typename T>
  TimeBasedBuckets<T>::TimeBasedBuckets(){
    /**/
  }

  // forward size query to underlying container
  template <typename T>
  size_t TimeBasedBuckets<T>::size() const{
    auto rv = _buckets.size();
    return rv;
  }

  // forward iterator access to underlying container
  template <typename T>
  TBBS_iterator<T> TimeBasedBuckets<T>::begin(){
    auto rv = _buckets.begin();
    return rv;
  }

  template <typename T>
  TBBS_const_iterator<T> TimeBasedBuckets<T>::begin() const{
    auto rv = _buckets.begin();
    return rv;
  }

  template <typename T>
  TBBS_iterator<T> TimeBasedBuckets<T>::end(){
    auto rv = _buckets.end();
    return rv;
  }

  template <typename T>
  TBBS_const_iterator<T> TimeBasedBuckets<T>::end() const{
    auto rv = _buckets.end();
    return rv;
  }

  // insertion
  template<typename T>
  void TimeBasedBuckets<T>::Insert(T item, double window){
    // if the item belongs in a preexisting bucket, put it there
    bool inserted = false;
    for (auto& bucket: (*this)){
      if (bucket.Accepts(item, window)){
        bucket.Append(item, window);
        inserted = true;
        return;
      }
    }

    // otherwise, this item defines a new bucket
    if (!inserted){
      TimeBasedBucket<T> bucket;
      _buckets.push_back(bucket);
      this->Insert(item, window);
    }
  }
}
