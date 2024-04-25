// Ed Callaghan
// Time-based bucketing, of e.g. StrawGasSteps and StrawDigis into
// groups of overlapping quasi-digitization windows
// February 2024

// stl
#include <vector>

namespace mu2e{
  // container-type which is aware of whether a candidate item should be
  // added, or not
  template<typename T>
  class TimeBasedBucket: public std::vector<T>{
    public:
      TimeBasedBucket(double);

      // whether a candidate item "fits" into this bucket's gross time window
      bool Accepts(T);
      // insertion --- assumes that appends are pre-sorted
      void Append(T);
    protected:
      /**/
    private:
      double window;
  };

  // a container of such containers, which ensures that every item considered
  // either lands in an existing bucket, or defines its own
  template <typename T>
  class TimeBasedBuckets: public std::vector< TimeBasedBucket<T> >{
    public:
      TimeBasedBuckets(double);
      // insertion
      void Insert(T);
    protected:
      /**/
    private:
      double window;
  };

  template <typename T>
  TimeBasedBucket<T>::TimeBasedBucket(double window): window(window){
    /**/
  }

  // whether a candidate item "fits" into this bucket's gross time window
  template <typename T>
  bool TimeBasedBucket<T>::Accepts(T candidate){
    // empty bucket rejects nothing
    if (this->size() < 1){
      return true;
    }

    // otherwise, check between first/last times, accounting for windowing
    double lower = this->front().time() - this->window;
    double upper = this->back().time() + this->window;
    double now = candidate.time();
    bool rv = (lower <= now) && (now <= upper);
    return rv;
  }

  // insertion --- assumes that appends are pre-sorted
  template<typename T>
  void TimeBasedBucket<T>::Append(T candidate){
    this->push_back(candidate);
  }

  template<typename T>
  TimeBasedBuckets<T>::TimeBasedBuckets(double window): window(window){
    /**/
  }

  // insertion
  template<typename T>
  void TimeBasedBuckets<T>::Insert(T item){
    // if the item belongs in a preexisting bucket, put it there
    bool inserted = false;
    for (auto& bucket: (*this)){
      if (bucket.Accepts(item)){
        bucket.Append(item);
        inserted = true;
        return;
      }
    }

    // otherwise, this item defines a new bucket
    if (!inserted){
      TimeBasedBucket<T> bucket(window);
      this->push_back(bucket);
      this->Insert(item);
    }
  }
}
