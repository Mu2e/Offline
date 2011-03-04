#ifndef maybe_ref_HH
#define maybe_ref_HH
//
// A "safe reference" class template to be used as a return type when it is possible
// that there is nothing to return.  
//
// The usual usage is expected to be something like:
//
//  cet::maybe_ref<T const> r(call some accessor);
//  if ( !t ) return;       // or continue, or break ....
//  T const& t(r.ref());
//   ... use t in your following code.
//
// $Id: maybe_ref.hh,v 1.1 2011/03/04 19:54:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/03/04 19:54:04 $
//
// Original author Rob Kutschke
//
//
// Notes:
// 1) For the usual Mu2e related reasons I want the accessor to be by const&
//    even though internally it holds a const*.  
//
// 2) This class will eventually move out of the Mu2e code base into the
//    cetlib library.  To minimize downstream maintenance, I have put it 
//    in the final namespace now.
//
// 3) Any pointer with a non-null value is presumed to be valid.
//
// 4) This class presumes that it is managing a non-owning pointer.  So it
//    will never call delete on the pointer and will always do shallow copies.
//
// 5) This class does not protect against accessing an object that
//    goes out of scope before this class goes out of scope, or
//    against accessing an object that moved in memory (such as in an
//    element in a std::vector after reallocation forced by a
//    push_back).  If a maybe_ref<T> is used within a method of a
//    framework module, any object obtained from the event or from a 
//    Service will satisfy this constraint.
//
// 6) An operator<() was not included on purpose. Sorting pointers is often
//    meaningless so we do not want to make it too easy.
//
// 7) The operator-> and operator* were not included because the desired usage is 
//    to immediately use the ref() method and use the resulting T const&.  We do not
//    want to tempt people to use this class like a pointer since it might get used
//    that way inside nested loops.
//
// 8) Do we want to make assignment private and unimplemented ( = delete in C++0X)?
//

#include <stdexcept>

namespace cet {

  template<typename T>
  class maybe_ref{

  public:

    typedef T value_type;

    // Default constructor makes an invalid object.
    maybe_ref():_ptr(0){}

    explicit maybe_ref( T* ptr):_ptr(ptr){}

    // Accept the compiler written d'tor, copy c'tor and assignment operator.

    bool isValid() const { return _ptr !=0; }

    // A synonym for isValid.
    operator bool() const { return _ptr !=0; }

    T& ref(){
      if ( !isValid() ) {
        throw std::logic_error(
          "Attempt to access invalid pointer using cet::maybe_ref ");
      }
      return *_ptr;
    }

    // For use when the optional<T> itself is const.
    T const& ref() const{
      if ( !isValid() ) {
        throw std::logic_error(
          "Attempt to access invalid pointer using cet::maybe_ref ");
      }
      return *_ptr;
    }

    void swap ( maybe_ref& rhs ){
      std::swap(_ptr, rhs._ptr);
    }

  private:

    // Non-owning pointer to the object.
    T* _ptr;

};

} // end namespace cet


#endif
