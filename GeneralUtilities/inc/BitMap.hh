
#ifndef GeneralUtilities_BitMap_hh
#define GeneralUtilities_BitMap_hh
//
// Template used to instantiate the bit map classes.
//
//
// The user must supply a detail class with the following requirements:
//
//   1) The detail class must contain an enum named bit_type.
//      The legal values of the enums are bit numbers that fit inside the mask_type, defined next.
//
//   2) The detail class must provide a typedef, named mask_type that specifices
//      the data type of the bitmap member datum.
//
//      If mask_type is char, then legal bit values are on the domain [0,7].
//      If mask_type is unsigned short, then legal bit values are on the domain [0,15].
//      And so on.  This has been tested for mask_type up to unsigned long (64 bits).
//
//   3) The detail class must contain three static functions:
//       static std::string const& typeName();
//       static std::map<std::string,mask_type> const& bitNames();
//       static mask_type bit_to_mask ( bit_type c );
//
//      The first function returns one string that holds the name of the type, usually the same name as the .hh file.
//      The template uses it decorating some printed output.
//
//      The second function returns an std::map that implements the cross reference between each mask_type value and
//      its string representation. Note that this is NOT a cross-reference between bit names and their
//      string representations.
//
//      The third function converts a bit_type to the corresponding mask_type.
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
// 3) By design there are no overloads of merge, clear, hasAllProperties or hasAnyProperties that take an argument
//    that is of built-in integral type.  This is needed for as a safety feature - see note 4 for an example.
//    The are also no overloads that take an argument of mask_type, to prevent people from hand constructing
//    illegal masks and using them.
//
// 4) Consider the following pieces of wrong code:
//
//       BitMap m;
//       m.merge(DETAIL::property1|DETAIL::property2);
//    or
//       m.hasAllProperties(DETAIL::property1|DETAIL::property2);
//
//    These will give erroneous results since it is not meaningful to OR objects of bit_type.  It is only legal to
//    OR objects of mask_type. These examples will not compile because result of the OR operation is a object of
//    type int and there are no overloads of merge, etc that take an int.
//
// 5) The code enforces the following notion of validity.  A value is valid if all of bits that are set in the
//    value correspond to bits that are defined by the bitNames method of the Detail class.  It is the responsibility
//    of the writer of the detail class to ensure that the enum declaration and the bitNames method are consistent.
//    It is possible to create an invalid value by casting an arbitrary int into a bit_type.
//

#include "GeneralUtilities/inc/toHex.hh"

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <regex>

namespace mu2e {

  template < class DETAIL > class BitMap : public DETAIL {

  public:

    typedef typename DETAIL::bit_type       bit_type;
    typedef typename DETAIL::mask_type      mask_type;
    typedef std::map<std::string,mask_type> map_type;

    explicit BitMap():_value(empty_value()){}

    // Explicit or not? See note 1
    BitMap( bit_type value): _value(DETAIL::bit_to_mask(value)) {}

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
      _value = DETAIL::bit_to_mask(c);
      return *this;
    }

    // Two merge methods; see notes 3 and 4.
    void merge( bit_type bitNumber){
      _value = static_cast<mask_type>( _value | DETAIL::bit_to_mask(bitNumber) );
    }

    void merge( BitMap arg){
      _value = static_cast<mask_type>( _value | arg._value);
    }

    BitMap& operator | (BitMap const& other) {
      merge(other);
      return *this;
    }

    void clear( bit_type bitNumber) {
      _value = static_cast<mask_type>(_value & ~DETAIL::bit_to_mask(bitNumber) );
    }

    void clear( BitMap arg) {
      _value = static_cast<mask_type>(_value & ~arg._value);
    }

    void reset(){
      _value = empty_value();
    }

    // By design there is no merge method taking an argument of built-in integral type.

    // Accessors.
    std::string hex() const { return toHex(static_cast<unsigned long>(_value)); }

    bool empty() const{
      return (_value == static_cast<mask_type>(0) );
    }

    // For convenience, forward the static function to a member function.
    bool isValid() const{
      return isValid(_value);
    }

    bool hasAllProperties( bit_type arg) const {
      return hasAllProperties(BitMap(arg));
    }

    bool hasAllProperties( BitMap arg) const {
      return (_value & arg._value ) == arg._value;
    }

    bool hasAnyProperty( BitMap arg) const {
      return (_value & arg._value ) != 0;
    }

    bool hasAnyProperty( bit_type arg) const {
      return hasAnyProperty(BitMap(arg));
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
        if ( hasAllPropertiesByMask(i->second) ) {
          if ( properties.size() > 0 ) properties += " ";
          properties += i->first;
        }
      }
      return properties;

    }

    // Implicit conversion versions of the accessors.
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

    bool lessByValue( BitMap arg ) const{
      return ( _value < arg._value );
    }

    // A value is invalid if any bits not defined in the detail class are set.
    static bool isValid( int arg ){
      static mask_type illegalBits = ~legalBits();
      return ( (arg & illegalBits) == 0);
    }

    static bool isValidOrThrow(int value){
      if ( !isValid(value) ){
        std::ostringstream os;
        os << DETAIL::typeName() << " invalid mask value : " << value;
        throw std::out_of_range( os.str() );
      }
      return true;
    }

    // The number of defined bits.
    static size_t size(){ return DETAIL::bitNames().size(); }

    // Access the translation map.
    static map_type const& bitNames() { return DETAIL::bitNames(); }

    // Print all values in the translation map.
    static void printAll( std::ostream& ost ){
      ost << "Defined mask values for " << DETAIL::typeName() << std::endl;
      for ( typename map_type::const_iterator i=bitNames().begin(), e=bitNames().end();
            i != e; ++i ){

        // Need the temporary in case mask_type is char; we want it to be treated as an unsigned.
        unsigned long tmp(i->second);

        ost << std::setw(10) << toHex(tmp) << " " << i->first << std::endl;
      }
    }

  private:

    mask_type _value;

    // Helper methods that expose mask_type.  Keep them private to force people to use
    // type safe(r) public interface.

    mask_type empty_value(){
      return static_cast<mask_type>(0);
    }

    bool hasAllPropertiesByMask( mask_type mask) const {
      return (_value & mask ) == mask;
    }

    mask_type findMaskByNameOrThrow( std::string const& name ){
      mask_type mask(0);
      bool invalid(false);
      const std::regex separator("[^\\s,:]+"); // begining or separator
      auto sbegin = std::sregex_iterator(name.begin(), name.end(), separator);
      auto send = std::sregex_iterator();
      for (std::sregex_iterator match = sbegin; match != send; ++match) {
	std::string subname = match->str();
	typename map_type::const_iterator j = bitNames().find(subname);
	if ( j == bitNames().end() ){
	  invalid = true;
	  break;
	}
	mask |= j->second;
      }
      if(invalid){
      // try to interpret as a (text) hex string; make sure to edit out spaces first!
	std::string cname(name);
	cname.erase(cname.begin(), std::find_if(cname.begin(), cname.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	if(cname.compare(0,2,"0x") == 0 || cname.compare(0,2,"0X") == 0){
	  mask = std::stoul(name,0,16);
	  if((isValid(mask)) && mask != 0){
	    invalid = false;
	  }
	}
	// if it's still invalid, throw
	if(invalid){
	  std::ostringstream os;
	  os << DETAIL::typeName() << " invalid mask name : " << name;
	  throw std::out_of_range( os.str() );
	}
      }
      return mask;
    }

    // Compute a mask in which all bits defined in the detail class are set.
    static mask_type legalBits(){

      mask_type tmp(0);
      map_type const& bitNames = DETAIL::bitNames();
      for ( typename map_type::const_iterator i=bitNames.begin(), e=bitNames.end();
            i != e; ++i ){
          tmp = tmp | i->second;
      }
      return tmp;
    }
  };

  template < class DETAIL >
  inline
  std::ostream& operator<<(std::ostream& ost,
                           const BitMap<DETAIL> & value ){

    std::string validity = value.isValid() ? "" :
      "Invalid value: undefined bits have been set: The known bits are: ";
    ost << "( "
        << value.hex() << " : "
        << validity
        << value.stringRep()
        << " )";
    return ost;
  }

  template < class DETAIL >
  inline
  bool lessByValue( BitMap<DETAIL> const& rhs, BitMap<DETAIL> const& lhs ){
    return ( rhs.lessByValue(lhs) );
  }

  template < class DETAIL >
  inline
  bool lessByStringRep( BitMap<DETAIL> const& rhs, BitMap<DETAIL> const& lhs ){
    return ( rhs.stringRep() < lhs.stringRep() );
  }

} // end namespace mu2e

#endif /* GeneralUtilities_BitMap_hh */
