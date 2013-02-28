#ifndef GeneralUtilities_BitMap_hh
#define GeneralUtilities_BitMap_hh
//
// Template used to instantiate the bit map classes.
//
//   $Id: BitMap.hh,v 1.2 2013/02/28 18:07:34 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2013/02/28 18:07:34 $
//
// The user must supply a detail class with the following requirements:
//
//   1) The detail class must contain an enum named bit_type.
//      The values of the enums are bit numbers on the domain [0,31].
//
//   2) The detail class must contain two static functions:
//       static std::string const& typeName();
//       static std::map<mu2e::BitMap_mask_type,std::string> const& bitNames();
//
//      The first function returns one string that holds the name of the type, usually the same name as the .hh file.
//      The template uses it decorating some printed output.
//
//      The second function returns an std::map that implements the cross reference between each mask_type value and
//      its string representation. Note that this is NOT a cross-reference between bit names and their string representations.
//
// Notes:
//
// 0) There is no c'tor from any primitive integral type; this is to enforce type safety so that you cannot cross-stitch
//    detail classes.
//
// 1) Do we want the c'tor from bit_type to be explicit or not? Consider the code:
//        Color v =  Color::red;
//    Do we want this to compile?  If yes, the the c'tor must be not explicit.  If we do not want it to compile
//    then it must be explicit.  A related question is illustrated by this code fragment:
//        Color v;
//        v = Color::Red;
//    If we want this to work, then we need and additional assignment operator; if do not want it to compile,
//    we do not want the additional assignment operator. In the present release, both examples work; we believe this
//    is safe because there are no types that can automatically convert to bit_type.
//
// 2) By design there is no operator<( BitMap const& ) since that would be ambiguous for comparison
//    by id or by string representation. Instead there are two free functions, lessByValue and lessByStringRep.
//
// 3) By design there are no overloads of merge or hasProperties that take an argument that is of built-in
//    integral type.  This is needed for as a safety feature - see note 4 for an example.
//    The are also no overloads that take an argument of mask_type, to prevent people from hand constructing
//    illegal masks and using them.
//
// 4) Consider the following:
//       BitMap m;
//       m.merge(DETAIL::property1|DETAIL::property2);
//    or
//       m.hasProperties(DETAIL::property1|DETAIL::property2);
//
//    These will give erroneous results since it is meaningful to OR objects of mask_type but not of bit_type.
//    These examples will not compile because result of the OR operation is a object of type int and there are
//    no overloads of merge or hasProperties that take an int.
//
// 5) The code enforces the following notion of validity.  A value is valid if all of bits that are set in the
//    value correspond to bits that are defined by the bitNames method of the Detail class.  It is the responsibility
//    of the writer of the detail class to ensure that the enum declaration and the bitNames method are consistent.
//    It is possible to create an invalid value by casting an arbitrary int into a bit_type.
//


#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace mu2e {

  // This type should be used by the Detail classes as the first template parameter of
  // the type returned by their bitNames method.
  typedef unsigned BitMap_mask_type;

  template < class DETAIL > class BitMap : public DETAIL {

  public:


    typedef typename DETAIL::bit_type  bit_type;
    typedef          BitMap_mask_type  mask_type;
    typedef std::map<mask_type,std::string> map_type;

    explicit BitMap():_value(empty_value()){}

    // Explicit or not? See note 1
    BitMap( bit_type value): _value(1<<value) {}

    // Constructor from a vector of bit names; must match the names used in the detail class.
    explicit BitMap( std::vector<std::string> const& names ):_value(empty_value()){

      for ( std::vector<std::string>::const_iterator i=names.begin();
            i != names.end(); ++i ){
        _value = _value | findMaskByNameOrThrow(*i);
      }

    }

    // Constructor from a string holding the name of bit.  Special case of the previous c'tor.
    explicit BitMap( std::string const& name ):_value(findMaskByNameOrThrow(name)){
    }

    // Accept compiler generated d'tor, copy c'tor and copy assignment.

    // Construct an object with all known bits set.
    static BitMap allBits(){
      BitMap b;
      b._value = legalBits();
      return b;
    }

    // Additional assignment operator. See note 1.
    BitMap& operator=( bit_type c){
      _value = 1<<c;
      return *this;
    }

    // Two merge methods; see notes 3 and 4.
    void merge( bit_type bitNumber){
      _value = static_cast<mask_type>( _value | (1<<bitNumber) );
    }

    void merge( BitMap arg){
      _value = static_cast<mask_type>( _value | arg._value);
    }

    void reset(){
      _value = empty_value();
    }

    // By design there is no merge method taking an argument of built-in integral type.

    // Accessors.
    mask_type value() const { return _value;}

    bool empty() const{
      return (_value == static_cast<mask_type>(0) );
    }

    // For convenience, forward the static function to a member function.
    bool isValid() const{
      return isValid(_value);
    }

    bool hasProperties( bit_type arg) const {
      return hasProperties(BitMap(arg));
    }

    bool hasProperties( BitMap arg) const {
      return (_value & arg._value ) == arg._value;
    }

    // Form a string representation of the bitmap.
    std::string stringRep() const {

      if ( empty() ){
        return "noBitsSet";
      }

      std::string properties;

      map_type const& bitNames = DETAIL::bitNames();
      for ( typename map_type::const_iterator i=bitNames.begin(), e=bitNames.end();
            i != e; ++i ){
        if ( hasPropertiesByMask(i->first) ) {
          if ( properties.size() > 0 ) properties += " ";
          properties += i->second;
        }
      }
      return properties;

    }

    // Implicit conversion versions of the accessors.
    operator mask_type ()const{
      return _value;
    }

    operator std::string() const{
      return stringRep();
    }

    // Comparisons. About operator< , see note 2.
    bool operator==( bit_type b) const{
      BitMap tmp(b);
      return ( _value == tmp._value );
    }

    bool operator==( BitMap arg) const{
      return ( _value == arg._value );
    }

    // A value is invalid if any bits not defined in the detail class are set.
    static bool isValid( int value ){
      static mask_type illegalBits = ~legalBits();
      return ( (value & illegalBits) == 0);
    }

    static bool isValidOrThrow(int value){
      if ( !isValid(value) ){
        std::ostringstream os;
        os << DETAIL::typeName() << " invalid mask value : " << value;
        throw std::out_of_range( os.str() );
      }
      return true;
    }

    // Tne number of defined bits.
    static size_t size(){ return DETAIL::bitNames().size(); }

    // Access the translation map.
    static map_type const& bitNames() { return DETAIL::bitNames(); }

    // Print all values in the translation map.
    static void printAll( std::ostream& ost ){
      ost << "Defined mask values for " << DETAIL::typeName() << std::endl;
      for ( typename map_type::const_iterator i=bitNames().begin(), e=bitNames().end();
            i != e; ++i ){
        ost << std::setw(10) << i->first << " " << i->second << std::endl;
      }
    }


  private:

    mask_type _value;

    // Helper methods that expose mask_type.  Keep them private to force people to use
    // type safe(r) public interface.

    mask_type empty_value(){
      return static_cast<mask_type>(0);
    }

    bool hasPropertiesByMask( mask_type mask) const {
      return (_value & mask ) == mask;
    }

    mask_type findMaskByNameOrThrow( std::string const& name ){

      typename map_type::const_iterator j = findMaskByName(name);
      if ( j == bitNames().end() ){
        std::ostringstream os;
        os << DETAIL::typeName() << " invalid mask name : " << name;
        throw std::out_of_range( os.str() );
      }
      return j->first;

    }

    typename map_type::const_iterator findMaskByName( std::string const& name ){

      for ( typename map_type::const_iterator j=bitNames().begin();
            j !=bitNames().end(); ++j ) {
        if ( j->second == name ) return j;
      }
      return bitNames().end();
    }

    // Compute a mask in which all bits defined in the detail class are set.
    static mask_type legalBits(){

      mask_type tmp(0);
      map_type const& bitNames = DETAIL::bitNames();
      for ( typename map_type::const_iterator i=bitNames.begin(), e=bitNames.end();
            i != e; ++i ){
          tmp = tmp | i->first;
      }

      return tmp;
    }

  };

  template < class DETAIL >
  inline
  std::ostream& operator<<(std::ostream& ost,
                           const BitMap<DETAIL> & value ){

    std::string validity = BitMap<DETAIL>::isValid(value) ? "" : "Invalid value: undefined bits have been set:  ";
    ost << "( "
        << value.value() << ": "
        << validity
        << value.stringRep()
        << " )";
    return ost;
  }

  template < class DETAIL >
  inline
  bool lessByValue( BitMap<DETAIL> const& rhs, BitMap<DETAIL> const& lhs ){
    return ( rhs.value() < lhs.value() );
  }

  template < class DETAIL >
  inline
  bool lessByStringRep( BitMap<DETAIL> const& rhs, BitMap<DETAIL> const& lhs ){
    return ( rhs.stringRep() < lhs.stringRep() );
  }

} // end namespace mu2e

#endif /* GeneralUtilities_BitMap_hh */
