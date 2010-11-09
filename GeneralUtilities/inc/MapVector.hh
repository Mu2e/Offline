#ifndef Generalutilities_MapVector_hh
#define Generalutilities_MapVector_hh
//
// An STL-like class template that looks and feels std::map<key_type,T>,
// with the exception that it is has a few extra modifier and accessor
// functions described below.  So far as I know, all accessor and 
// modifier methods of std::map work; I have not implemented the
// constructors that allow user specified comparator and allocator objects.
//
// $Id: MapVector.hh,v 1.8 2010/11/09 20:53:58 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/09 20:53:58 $
//
//   Original author Rob Kutschke
//
//  The additional methods are:
//  VALUE&       findOrThrow( key_type key );
//      - If the key exists, it returns a reference to the value.
//      - If the key does not exist it throws.
//
//  VALUE const& findOrThrow( key_type key ) const;
//      - same as the previous method except that it may run on a 
//        const MapVector and returns a const reference.
//
//  VALUE*       findOrNull( key_type key );
//      - If the key exists, it returns a non-owning pointer to the value.
//      - If the key does not exist it returns a null pointer.
//
//  VALUE const* findOrNull( key_type key ) const;
//      - same as the previous method except that it may run on a 
//        const MapVector and returns a const pointer.
//
//  VALUE const& operator[]( key_type key ) const;
//      - Alternate syntax for the previous method.
//      - This is the original use case for this class.
//
//  VALUE const& at( key_type key ) const;
//      - Alternate syntax for the previous method.
//      - Needed since existing user code uses it
//
//   bool has( key_type key ) const;
//      - Return true if the key exists in the map.
//      - Otherwise return false.
//
//  iterator insertOrThrow ( const value_type& value );
//      - If the key does not exist in the map, insert the pair and return
//        an iterator to the pair.
//      - If the key already exists in the map, throw.
//
// Usage notes:
//  The ONLY way to write a loop over a MapVector is:
//      MapVector<T> const& mvt = ....;
//      for ( MapVector<T>::const_iterator i=mvt->begin();
//             i!=mvt->end(); ++i ){
//           T const& t = i->second;
//      }
//
//  If you write a loop like:
//      MapVector<T> const& col = ....;
//      for ( size_t i=0; i<col.end(); ++i){
//        T const& t = col[i].
//      }
//  This will give a compiler error because there is no operator[] that takes
//  an int as an argument, only an operator[] that takes a 
//  MapVector<T>::key_type as an argument.
//
// Other Notes:
// This class is designed to address the following use case:
//
// 1) I have some objects of type T, each of which has an identifier,
//    of type key_type, where key_type is an unsigned integral type.
//
// 2) The key_type is not dense; for example, it might be legal to
//    have three objects with identifiers with of { 1, 99, 4567 }.
//
// 3) I would like to make a collection of these with a operator[]
//    accessor method that works on a const collection, and which
//    throws if the identifier does not exist.
//
// 4) I cannot (efficiently) use std::vector<T> since the identifiers are not dense.
//    I cannot use std::map<key_type,T> since it fails point 3).
//
// 5) So I have invented this class that provides a const operator []
//    and other wise forwards to std::map.  In the process I added a few
//    other methods.
//
// 6) In the future we can replace the underlying container with 
//    a different implementation that is better optimized for lookup.
//
// Todo:
//  1) Add c'tor allowing use specified Compare and Allocator??
//  2) Are there better names for findOrThrow and insertOrThrow
//  3) Which methods should I comment out since they are not likely to be in the final MapVector?
//  4) Any other methods that I want?
//  5) Do I want to hide sstream to avoid it being in the header?
//
#include <map>
#include <stdexcept>
#include <sstream>

// Mu2e includes
#include "GeneralUtilities/inc/MapVectorKey.hh"

template<typename VALUE>
class MapVector{

public:

  typedef MapVectorKey                      key_type;
  typedef VALUE                          mapped_type;
  typedef std::pair<key_type,mapped_type> value_type;
  typedef std::map<key_type,VALUE>           MapType;

  typedef typename MapType::size_type             size_type;
  typedef typename MapType::difference_type difference_type;
  
  typedef typename MapType::iterator                             iterator;
  typedef typename MapType::const_iterator                 const_iterator;
  typedef typename MapType::reverse_iterator             reverse_iterator;
  typedef typename MapType::const_reverse_iterator const_reverse_iterator;

  typedef typename MapType::key_compare     key_compare;
  typedef typename MapType::value_compare value_compare;

  typedef typename MapType::allocator_type                   allocator_type;
  typedef typename MapType::allocator_type::pointer                 pointer;
  typedef typename MapType::allocator_type::const_pointer     const_pointer;
  typedef typename MapType::allocator_type::reference             reference;
  typedef typename MapType::allocator_type::const_reference const_reference;

  // Constructors
  MapVector():_map(){}

  template <class InputIterator>
  MapVector ( InputIterator first, InputIterator last):_map(first,last){}

  // Compiler generated d'tor, copy c'tor and assignment operator are OK.

  // The additional functions are first:

  // Return a reference to the value corresponding to an existing key,
  // throw if the key does not exist.
  VALUE& findOrThrow( key_type key ){
    iterator i = _map.find(key);

    if ( i == _map.end() ){
      std::ostringstream out;
      out << "No such key: " << key;
      throw std::out_of_range( out.str() );
    }
    return i->second;
  }

  // The const version of the above.
  VALUE const& findOrThrow( key_type key ) const{
    const_iterator i = _map.find(key);

    if ( i == _map.end() ){
      std::ostringstream out;
      out << "No such key: " << key;
      throw std::out_of_range( out.str() );
    }
    return i->second;
  }

  // Return a non-owning pointer to the value corresponding to an existing key,
  // Return a null pointer if the key does not exist.
  VALUE* findOrNull( key_type key ){
    iterator i = _map.find(key);

    if ( i == _map.end() ){
      return 0;
    }
    return &(i->second);
  }

  // Const version of the above.
  VALUE const * findOrNull( key_type key ) const{
    const_iterator i = _map.find(key);

    if ( i == _map.end() ){
      return 0;
    }
    return &(i->second);
  }

  // Two handy synonyms for the const findOrThrow;
  // These methods are the original reason for this class.
  VALUE const& operator[] ( key_type key ) const{
    return findOrThrow(key);
  }

  VALUE const& at( key_type key ) const{
    return findOrThrow(key);
  }

  // Check to see if the key exists in the map.
  bool has( key_type key ) const{
    const_iterator i = _map.find(key);
    return ( i != _map.end() );
  }

  // Insert the requested pair or throw if the key already exists.
  iterator insertOrThrow ( const value_type& value ){
    const_iterator i = _map.find(value.first);
    
    if ( i != _map.end() ){
      std::ostringstream out;
      out << "Found key: " << value.first << " when no key was expected.";
      throw std::out_of_range( out.str() );
    }
    return _map.insert(value).first;
  }

  // Everything below here just forwards to std::map.

  // Modifiers:
  VALUE& operator[] ( key_type key ){
    return _map[key];
  }

  std::pair<iterator,bool> insert ( const value_type& value ){
    return _map.insert(value);
  }

  template <class InputIterator>
  iterator insert ( InputIterator position, const value_type& value ){
    return _map( position, value);
  }

  void erase ( iterator position ){
    _map.erase(position);
  }

  size_type erase ( key_type key ){
    return _map.erase(key);
  }

  void erase ( iterator first, iterator last ){
    _map.erase(first,last);
  }

  void swap ( MapVector<VALUE>& mv ){
    _map.swap(mv._map);
  }

  void clear(){
    _map.clear();
  }

  // Capacity
  bool      empty()    const { return _map.empty();    }
  size_type size()     const { return _map.size();     }
  size_type max_size() const { return _map.max_size(); }

  // Iterators
  iterator begin (){
    return _map.begin();
  };

  const_iterator begin () const{
    return _map.begin();
  }

  iterator end (){
    return _map.end();
  };

  const_iterator end () const{
    return _map.end();
  }

  reverse_iterator rbegin (){
    return _map.rbegin();
  };

  const_reverse_iterator rbegin () const{
    return _map.rbegin();
  }

  reverse_iterator rend (){
    return _map.rend();
  };

  const_reverse_iterator rend () const{
    return _map.rend();
  }

  // Operations
  iterator find( key_type key) {
    return _map.find(key);
  }

  const_iterator find( key_type key) const{
    return _map.find(key);
  }

  size_type count ( key_type key ) const{
    return _map.count(key);
  }

  iterator lower_bound ( key_type key ){
    return _map.lower_bound(key);
  }

  const_iterator lower_bound ( key_type key ) const{
    return _map.lower_bound(key);
  }

  iterator upper_bound ( key_type key ){
    return _map.upper_bound(key);
  }

  const_iterator upper_bound ( key_type key ) const{
    return _map.upper_bound(key);
  }

  std::pair<iterator,iterator> equal_range ( key_type key ){
    return _map.equal_range(key);
  }

  std::pair<const_iterator,const_iterator> equal_range ( key_type key ) const{
    return _map.equal_range(key);
  }

  key_compare key_comp ( ) const{
    return _map.key_comp();
  }

  value_compare value_comp ( ) const{
    return _map.value_comp();
  }

  allocator_type get_allocator() const{
    return _map.get_allocator();
  }

private:

  std::map<key_type,VALUE> _map;

};

#endif
