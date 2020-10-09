#ifndef Sandbox_KalRepCollection_hh
#define Sandbox_KalRepCollection_hh
//
// A draft of a DataProduct class to take ownership of KalRep objects, provide
// access to them and manage their lifetimes.  This class is designed to be
// a transient only data product.
//
//
// Original author Rob Kutschke
//
// Notes.
// 1) Consider the following fragment of a producer method:
//    unique_ptr<T> t(new T);
//    // Fill t.
//    event.put(std::move(t));
//    In side event.put the following happens:
//       a) Default construct a T owned by the event. Call this object t0.
//       b) swap the contents of t and t0 using T::swap member function.
//       c) 
// 

// Forward reference
class KalRep;

namespace mu2e {

  class KalRepCollection{

  public:

    KalRepCollection():tracks_(){}

    // Caller transfers ownership of the KalRep objects to us.
    explicit KalRepCollection( std::vector<KalRep*> v ):tracks_(v){}

    ~KalRepCollection(){
      for( std::vector<KalRep*>::iterator i=tracks_.begin();
           i!=tracks_.end(); ++i ){
        KalRep* t = *i;
        if ( t ) delete t;
      }
      // tracks_ immediately goes out of scope:
    }

    // Caller transfers ownership of the KalRep object to us.
    void push_back( KalRep* t){
      tracks_.push_back(t);
    }

    void pop_back( ){
      KalRep* t = tracks_.back();
      if ( t ) delete t;
      tracks_.pop_back(t);
    }

    // Needed for event.put().
    void swap( KalRepCollection& rhs){
      std::swap( this->tracks_, rhs.tracks_);
    }

    // Accessors: this class retains ownership.
    KalRep const& operator[](size_t i) const { return *tracks_.at(i) }
    KalRep&       operator[](size_t i)       { return *tracks_.at(i) }
    KalRep const& at        (size_t i) const { return *tracks_.at(i) }
    KalRep&       at        (size_t i)       { return *tracks_.at(i) }

    // Or we could provide these accessors instead.
    //KalRep const* operator[](size_t i) const { return tracks_.at(i) }
    //KalRep*       operator[](size_t i)       { return tracks_.at(i) }
    //KalRep const* at        (size_t i) const { return tracks_.at(i) }
    //KalRep*       at        (size_t i)       { return tracks_.at(i) }

    // Also forward the front and back functions?
    // Do we provide iterators?

  private:

    // Not copyable or assignable.
    KalRepCollection( KalRepCollection const& );
    KalRepCollection& operator=( KalRepCollection const& );

    // Owning pointers to the KalRep objects.
    std::vector<KalRep*> tracks_;

  };

} // namespace mu2e

#endif /* Sandbox_KalRepCollection_hh */
