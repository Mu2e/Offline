#ifndef GeneralUtilities_BitMap_hh
#define GeneralUtilities_BitMap_hh
//
// Template used to instantiate the bit map classes.
//
//   $Id: BitMap.hh,v 1.1 2013/02/22 00:36:06 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2013/02/22 00:36:06 $
//
// The user must supply a detail class with the following requirements:
//
//   1) The class must contain an enum named mask_type.
//   2) The class must contain two static functions:
//       static std::string const& typeName();
//       static std::map<mask_type,std::string> const& bitNames();
//
//      The first function returns one string that holds the name of the type, usually the same name as the .hh file.
//      The template uses it decorating some printed output.
//
//      The second function returns an std::map that implements the cross reference between each mask_type value and
//      its string representation.
//
// Notes:
//
// 0) Do we want to remove c'tor from an int?  It breaks the type safety and the validity check will not catch it.
//    The downside is that there are valid reasons for to have it. Mabye a static conversion function makes more sense
//    since it needs to be used more consciously?
//
// 1) Do we want the c'tor from mask_type to be explicit or not? Consider the code:
//        TrackBits v =  TrackBits::passedStep2;
//    Do we want this to compile?  If yes, the the c'tor must be not explicit.  If we do not want it to compile
//    then it must be explicit.  A related question is illustrated by this code fragment:
//        TrackBits v;
//        v = TrackBits::passedStep2;
//    If we want this to work, then we need to additional assignment operator; if do not want it to compile,
//    we do not want the additional assignment operator.
//
// 2) By design there is no operator<( BitMap const& ) since that would be ambiguous for comparison
//    by id or by string representation. Instead there are two free functions, lessByValue and lessByStringRep.
//
// 3) By design there are no overloads of merge or hasProperties that take an argument that is of built-in
//    integral type.  This is needed to ensure type saftey: if I have an object of one bitmap type and use
//    enum values from a different bitmask type, I will get a compiler error.
//
//    The class template does provide a c'tor from an int; but it is slow because it checks that only known bits are set.
//
// 4) Consider the following:
//       BitMap m;
//       m.merge(DETAIL::property1|DETAIL::property2);
//    or
//       m.hasProperties(DETAIL::property1|DETAIL::property2);
//
//    In both cases, the expression in the argument returns a integral type; but, by note 3), there are no
//    overloads of merge or hasProperties that takes an integral type.  The solution is to provide an overload
//    of operator| that will evaluate the expression and return a value of type mask_type.
//
// 5) The code enforces the following notion of validity.  A value is valid if all of bits that are set in the
//    value correspond to bits that are defined by the bitNames method of the Detail class.  It is the responsibility
//    of the writer of the detail class to ensure that the enum declaration and the bitNames method are consistent.
//    It is possible to create an invalid value by casting an arbitrary int into a mask_type.
//
// 6) Do we want a second template parameter for the type of the value?  So what we can use short int or char when
//    that would cover the anticipated need for bits?  Maybe this becomes moot in C+11 if we use enum classes.
//


#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace mu2e {

  // Helper function.  See note 4.
  template <class T>
  T operator|(T a, T b)
  {
    return static_cast<T>(static_cast<int>(a) | static_cast<int>(b));
  }

  template < class DETAIL > class BitMap : public DETAIL {
  public:

    typedef typename DETAIL::mask_type      mask_type;
    typedef std::map<mask_type,std::string> map_type;

    explicit BitMap():_value(static_cast<mask_type>(0)){}

    // Explicit or not? See note 1
    BitMap( mask_type value): _value(value) {}

    // Constructor from a vector of bit names; must match the names used in the detail class.
    explicit BitMap( std::vector<std::string> const& names ):_value(static_cast<mask_type>(0)){

      for ( std::vector<std::string>::const_iterator i=names.begin();
            i != names.end(); ++i ){
        _value = _value | findMaskByNameOrThrow(*i);
      }

    }

    // Constructor from a string holding the name of bit.  Special case of the previous c'tor.
    explicit BitMap( std::string const& name ):_value(findMaskByNameOrThrow(name)){
    }

    // Accept compiler generated d'tor, copy c'tor and copy assignment.

    // Additional assignment operator. See note 1.
    BitMap& operator=( mask_type c){
      _value = c;
      return *this;
    }

    // Two merge methods; see notes 3 and 4.
    void merge( mask_type mask){
      _value = static_cast<mask_type>(_value | mask);
    }

    void merge( BitMap mask){
      merge( mask.value() );
    }

    void reset(){
      _value = static_cast<mask_type>(0);
    }

    // By design there is no merge method taking an argument of built-in integral type.

    // Accessors.
    mask_type value() const { return _value;}

    bool empty() const{
      return (_value == static_cast<mask_type>(0));
    }

    // For convenience, forward the static function to a member function.
    bool isValid() const{
      return isValid(_value);
    }

    // Two hasProperties methods; see notes 3 and 4.
    bool hasProperties( mask_type mask) const {
      return (_value & mask ) == mask;
    }

    bool hasProperties( BitMap mask) const {
      return hasProperties( mask.value());
    }

    // Form a string representation of the bitmap.
    std::string stringRep() const {

      if ( _value == static_cast<mask_type>(0) ){
        return "noBitsSet";
      }

      std::string properties;

      map_type const& bitNames = DETAIL::bitNames();
      for ( typename map_type::const_iterator i=bitNames.begin(), e=bitNames.end();
            i != e; ++i ){
        if ( hasProperties(i->first) ) {
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

    // Comparisons.  See note 2.
    bool operator==( BitMap g) const{
      return ( _value == g._value );
    }

    bool operator==( mask_type g) const{
      return ( _value == g );
    }

    // A value is invalid if any bits not defined in the detail class are set.
    static bool isValid( int value ){
      static int illegalBits = ~legalBits();
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

    // Helper method for c'tor
    mask_type findMaskByNameOrThrow( std::string const& name ){

      typename map_type::const_iterator j = findMaskByName(name);
      if ( j == bitNames().end() ){
        std::ostringstream os;
        os << DETAIL::typeName() << " invalid mask name : " << name;
        throw std::out_of_range( os.str() );
      }
      return j->first;

    }

    // Helper method for c'tor
    typename map_type::const_iterator findMaskByName( std::string const& name ){

      for ( typename map_type::const_iterator j=bitNames().begin();
            j !=bitNames().end(); ++j ) {
        if ( j->second == name ) return j;
      }
      return bitNames().end();
    }

    // Compute a mask in which all bits defined in the detail class are set.
    static int legalBits(){

      int tmp(0);
      map_type const& bitNames = DETAIL::bitNames();
      for ( typename map_type::const_iterator i=bitNames.begin(), e=bitNames.end();
            i != e; ++i ){
          tmp = tmp | i->first;
      }

      return tmp;
    }

  private:

    mask_type _value;

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
