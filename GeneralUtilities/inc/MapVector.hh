#ifndef MapVector_HH
#define MapVector_HH
//
// An STL-like class template that looks and feels std::map<size_t,T>,
// with the exception that it is as a few extra modifier and accessor
// functions described below.  So far as I know, all methods of std::map
// work.
//
// $Id: MapVector.hh,v 1.1 2010/11/05 15:19:21 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/05 15:19:21 $
//
//   Original author Rob Kutschke
//
//  The additional methods are:
//  VALUE&       findOrThrow   ( key_type key );
//      - If the key exists, it returns a reference to the value.
//      - If the key does not exist it throws.
//  VALUE const& findOrThrow   ( key_type key ) const;
//      - same as the previous method except that it may run on a 
//        const MapVector and returns a const reference.
//  VALUE const& operator[]    ( key_type key ) const;
//      - Alternate syntax for the previous method.
//      - This is the original use case for this class.
//  iterator     insertOrThrow ( const value_type& value );
//      - If the does not exist in the map, insert the pair and return
//        an iterator to the pair
//      - If the key already exists in the map, throw.
//
// Notes:
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
//
// Todo:
//  1) Fix c'tor syntax with Allocator and Compare.
//  2) do I really want insertOrThrow?
//  3) Are there better names for findOrThrow and insertOrThrow
//  4) Is size_t the right key_type (root persistence?) ?
//  5) Which methods should I comment out since they are not likely to be in the final MapVector?
//  6) Any other methods that I want?
//  7) Do I want to hide sstream to avoid it being in the header?
//
#include <map>
#include <stdexcept>
#include <sstream>

template<typename VALUE>
class MapVector{

public:

  typedef std::size_t                       key_type;
  typedef VALUE                          mapped_type;
  typedef std::pair<key_type,mapped_type> value_type;
  typedef std::map<key_type,VALUE>           MapType;

  typedef typename MapType::size_type             size_type;
  typedef typename MapType::difference_type difference_type;

  typedef typename MapType::iterator                             iterator;
  typedef typename MapType::const_iterator                 const_iterator;
  typedef typename MapType::reverse_iterator             reverse_iterator;
  typedef typename MapType::const_reverse_iterator const_reverse_iterator;

  typedef typename MapType::key_compare key_compare;
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

  // Not able to get these to work:
  //explicit MapVector ( const Compare& comp = Compare(),
  // const Allocator& = Allocator() ):_map(){}
  //MapVector ( InputIterator first, InputIterator last,
  //    const Compare& comp = Compare(), const Allocator& = Allocator() );

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

  // A handy synonym for the const findOrThrow;
  // This method is the original reason for this class.
  VALUE const& operator[] ( key_type key ) const{
    return findOrThrow(key);
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

  iterator insert ( iterator position, const value_type& value ){
    return _map( position, value);
  }

  void erase ( iterator position ){
    _map.erase(position);
  }

  size_type erase ( key_type key ){
    _map.erase(key);
  }

  void erase ( iterator first, iterator last ){
    _map.erase(first,last);
  }

  void swap ( MapVector& mv ){
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
    return _map.keycomp();
  }

  value_compare value_comp ( ) const{
    return _map.valuecomp();
  }

  allocator_type get_allocator() const{
    return _map.get_allocator();
  }

private:

  std::map<key_type,VALUE> _map;

};

#endif
