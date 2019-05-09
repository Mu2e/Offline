#ifndef GeneralUtilities_OwningPointerCollection_hh
#define GeneralUtilities_OwningPointerCollection_hh
//
// A class template to take ownership of a collection of bare pointers to
// objects, to provide access to those objects and to delete them when the
// container object goes out of scope.
//
// This is designed to allow complex objects made on the heap to be used
// as transient-only data products.
//
// The original use is for BaBar tracks.
//
// $Id: OwningPointerCollection.hh,v 1.9 2014/04/18 16:39:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/04/18 16:39:30 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) I believe it would be wise to make the underlying container:
//      std::vector<boost::right_kind_of_smart_pointer<T>>.
//    scoped_ptr will not work since they are not copyable ( needed for vector ).
//    unique_ptr seems mostly concerned with threads.
//    I tried shared_ptr but I could not make it play nice with genreflex.
//

#include <vector>
#include <memory>

#include "canvas/Persistency/Common/detail/maybeCastObj.h"

namespace mu2e {

  template<typename T>
  class OwningPointerCollection{

  public:

    typedef typename std::vector<const T*>::const_iterator const_iterator;
    typedef T const* value_type;

    OwningPointerCollection():v_(){}

    // Caller transfers ownership of the pointees to us.
    explicit OwningPointerCollection( std::vector<T*>& v ):
      v_(v.begin(), v.end() ){
    }

    // Caller transfers ownership of the pointees to us.
    explicit OwningPointerCollection( std::vector<value_type>& v ):
      v_(v){
    }

    // We own the pointees so delete them when our destructor is called.
    ~OwningPointerCollection(){
      for( typename std::vector<value_type>::iterator i=v_.begin();
           i!=v_.end(); ++i ){
        delete *i;
      }
    }

#ifndef __GCCXML__

    // GCCXML does not know about move syntax - so hide it.

    OwningPointerCollection( OwningPointerCollection && rhs ):
      v_(std::move(rhs.v_)){
      rhs.v_.clear();
    }

    OwningPointerCollection& operator=( OwningPointerCollection && rhs ){
      v_ = std::move(rhs.v_);
      rhs.v_.clear();
      return *this;
    }

#endif /* GCCXML */

    // Caller transfers ownership of the pointee to us.
    void push_back( T* t){
      v_.push_back(t);
    }

    void push_back( value_type t){
      v_.push_back(t);
    }

    typename std::vector<value_type>::const_iterator begin() const{
      return v_.begin();
    }

    typename std::vector<value_type>::const_iterator end() const{
      return v_.end();
    }

    typename std::vector<value_type>::const_iterator cbegin() const{
      return v_.cbegin();
    }

    typename std::vector<value_type>::const_iterator cend() const{
      return v_.cend();
    }

    // Possibly needed by producer modules?
    /*
    void pop_back( ){
      delete v_.back();
      v_.pop_back();
    }
    */

    // Needed for event.put().
    void swap( OwningPointerCollection& rhs){
      std::swap( this->v_, rhs.v_);
    }

    // Accessors: this container retains ownership of the pointees.
    size_t     size()                    const { return   v_.size(); }
    T const&   operator[](std::size_t i) const { return  *v_.at(i);  }
    T const&   at        (std::size_t i) const { return  *v_.at(i);  }
    value_type get       (std::size_t i) const { return   v_.at(i);  }
    value_type operator()(std::size_t i) const { return   v_.at(i); }

    // const access to the underlying container.
    std::vector<value_type> const& getAll(){ return v_; }

  private:

    // Not copy-copyable or copy-assignable; this is needed to ensure exactly one delete.
    // GCCXML does not know about =delete so leave these private and unimplemented.
    OwningPointerCollection( OwningPointerCollection const& );
    OwningPointerCollection& operator=( OwningPointerCollection const& );

    // Owning pointers to the objects.
    std::vector<value_type> v_;

  };

} // namespace mu2e



// Various template specializations needed to make an art::Ptr<T> into an OwningPointerCollection<T>
// work.
//
// ItemGetter          - return a bare pointer to an requested by giving its index.
// has_setPtr          - do specializations exists for setPtr and getElementAddresses
// setPtr              - return a bare pointer to an requested by giving its index.
// getElementAddresses - return a vector of bare pointers to a vector of elements requested by index.
//
#ifndef __GCCXML__

#include "canvas/Persistency/Common/Ptr.h"

namespace art {
  namespace detail {
    template <typename T>
    class ItemGetter<T, mu2e::OwningPointerCollection<T> >;
  }
  template <class T>
  struct has_setPtr<mu2e::OwningPointerCollection<T> >;
}

namespace mu2e {
  template <class T>
  void
  setPtr(OwningPointerCollection<T> const & coll,
         const std::type_info & iToType,
         unsigned long iIndex,
         void const *& oPtr);

  template <typename T>
  void
  getElementAddresses(OwningPointerCollection<T> const & obj,
                      const std::type_info & iToType,
                      const std::vector<unsigned long>& iIndices,
                      std::vector<void const *>& oPtr);

}

template <typename T>
class art::detail::ItemGetter<T, mu2e::OwningPointerCollection<T> > {
public:
  T const * operator()(mu2e::OwningPointerCollection<T> const * product,
                       typename art::Ptr<T>::key_type iKey) const;
};

template <typename T>
inline
T const *
art::detail::ItemGetter<T, mu2e::OwningPointerCollection<T> >::
operator()(mu2e::OwningPointerCollection<T> const * product,
           typename art::Ptr<T>::key_type iKey) const
{
  assert(product != 0);
  std::size_t i(iKey);
  return product->get(i);
}

namespace art {
  template <class T>
    struct has_setPtr<mu2e::OwningPointerCollection<T> >
  {
    static bool const value = true;
  };
}

namespace mu2e {
  template <class T>
  void
  setPtr(OwningPointerCollection<T> const & coll,
         const std::type_info & iToType,
         unsigned long iIndex,
         void const *& oPtr)
  {
    oPtr = art::detail::maybeCastObj(coll.get(iIndex),iToType);
  }

  template <typename T>
  void
  getElementAddresses(OwningPointerCollection<T> const & obj,
                      const std::type_info & iToType,
                      const std::vector<unsigned long>& iIndices,
                      std::vector<void const *>& oPtr)
  {
    oPtr.reserve(iIndices.size());
    for ( auto i : iIndices ){
      oPtr.push_back(art::detail::maybeCastObj(obj.get(i),iToType));
    }

  }

}

#endif // __GCCXML__

#endif /* GeneralUtilities_OwningPointerCollection_hh */
