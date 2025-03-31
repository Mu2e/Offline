#ifndef GeneralUtilities_EnumToStringSparse_hh
#define GeneralUtilities_EnumToStringSparse_hh
//
// Template used to instantiate enum-matched-to-string classes.
//
//
// The user must supply a Detail class with the following requirements:
//
//   1) The class must contain an enum named enum_type.
//   2) The enum_type must contain a value named "unknown".  See note 1, below.
//   3) The class must contain two static functions:
//       static std::string const& typeName();
//       static std::map<enum_type,std::string> const& names();
//
//      The first function returns a string that holds the name of the type.
//      The template uses this function to decorate some printed output.
//      The suggested practice is as follows:  if the type ColorId is defined as:
//        typedef EnumToStringSparse<ColorIdDetail> ColorId;
//      then the method ColorIdDetail::typename should return the string "ColorId".
//
//      The second function returns a std::map that implements the cross reference between the enum_type and
//      the string representations of the enum value.
//
// Notes:
//
// 1) The ROOT IO system and several useful STL container types require that the class have a
//    default c'tor.  Therefore we require the enum_type to have a value named Detail::unknown .
//    Default constructed objects have this value.
//
// 2) By design there is no operator<( EnumToStringSparse const& ) since that would be ambiguous for comparison
//    by id or by name. Instead there are two free functions, lessById and lessByName.
//
// 3) There are two notions of validity.
//      a) The value is defined by the enum type and the value Detail::unknown is considered valid.
//      b) The value is defined by the enum type but the value Detail::unknown is NOT considered valid.
//    Several functions have a second argument that permits the user to choose between these two notions.
//    The default behaviour is the value Detail::unknown is NOT considered valid.
//

#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace mu2e {

  template < class Detail > class EnumToStringSparse : public Detail {
  public:

    typedef typename Detail::enum_type                       enum_type;
    typedef typename Detail::enum_type                       type;
    typedef std::map<typename Detail::enum_type,std::string> map_type;

    // See note 1.
    explicit EnumToStringSparse():_id(Detail::unknown){}

    explicit EnumToStringSparse( typename Detail::enum_type id): _id(id) {}

    explicit EnumToStringSparse( int id, bool unknownIsError=true ):
      _id(static_cast<typename Detail::enum_type>(id)){
      isValidOrThrow(id,unknownIsError);
    }

    explicit EnumToStringSparse ( std::string const& name,
                            bool throwIfUnknown=true,
                            bool throwIfUndefined=true):
      _id(findByName(name,throwIfUnknown,throwIfUndefined).id()){
    }

    // Accept compiler generated d'tor, copy c'tor and copy assignment.

    // One more assignment operator.
    EnumToStringSparse& operator=(typename Detail::enum_type c){
      _id = c;
      return *this;
    }

    // Accessors.
    typename Detail::enum_type id() const { return _id;}

    std::string const& name() const {

      // We checked for the validity of _id upon construction.  No need to check now.
      typename map_type::const_iterator value = Detail::names().find(_id);
      return value->second;
    }

    // Implicit conversion versions of the accessors.
    operator typename Detail::enum_type ()const{
      return _id;
    }

    operator std::string const& () const{
      return name();
    }

    // Comparisons.  See note 2.
    bool operator==(const EnumToStringSparse g) const{
      return ( _id == g._id );
    }

    bool operator==(const typename Detail::enum_type g) const{
      return ( _id == g );
    }

    // See note 3.
    static bool isValid( int id, bool unknownIsError=true ){

      typename Detail::enum_type cid = static_cast<typename Detail::enum_type>(id);
      typename map_type::const_iterator it = Detail::names().find(cid);

      if ( it == names().end() ) return false;

      return unknownIsError ? (cid !=Detail::unknown) : true;

    }

    static bool isValidOrThrow(int id, bool unknownIsError=true ){
      if ( !isValid(id,unknownIsError) ){
        std::ostringstream os;
        os << Detail::typeName() << " invalid enum value : " << id;
        throw std::out_of_range( os.str() );
        return false;
      }
      return true;
    }

    // Find an enum-value by its name; see note 3.
    static EnumToStringSparse findByName( std::string const& name,
                                    bool throwIfUnknown=true,
                                    bool throwIfUndefined=true ){

      map_type const& n = Detail::names();
      for ( typename map_type::const_iterator i=n.begin(), e=n.end();
            i != e; ++i ){
        if ( i->second == name ){
          if ( i->first == Detail::unknown && throwIfUnknown ){
            std::ostringstream os;
            os << Detail::typeName() << "::unknown is not allowed at this time";
            throw std::out_of_range( os.str() );
          }
          return EnumToStringSparse(i->first);
        }
      }

      // Did not find it;
      if ( throwIfUndefined ){
        std::ostringstream os;
        os << Detail::typeName() << " invalid enum name : " << name;
        throw std::out_of_range( os.str() );
      }

      return EnumToStringSparse();
    }

    // Tne number of names, including "unknown".
    static std::size_t size(){ return Detail::names().size(); }

    // Access the translation map.
    static map_type const& names() { return Detail::names(); }

    static void printAll( std::ostream& ost ){
      constexpr unsigned fieldWidth{10};
      ost << "Defined enum values for " << Detail::typeName() << std::endl;
      for ( typename map_type::const_iterator i=names().begin(), e=names().end();
            i != e; ++i ){
        ost << std::setw(fieldWidth) << i->first << " " << i->second << std::endl;
      }
    }

  private:

    typename Detail::enum_type _id;

  };

  template < class Detail >
  std::ostream& operator<<(std::ostream& ost,
                           const EnumToStringSparse<Detail> & id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

  template < class Detail >
  bool lessById( EnumToStringSparse<Detail> const& rhs, EnumToStringSparse<Detail> const& lhs ){
    return ( rhs.id() < lhs.id() );
  }

  template < class Detail >
  bool lessByName( EnumToStringSparse<Detail> const& rhs, EnumToStringSparse<Detail> const& lhs ){
    return ( rhs.name() < lhs.name() );
  }

}

#endif /* GeneralUtilities_EnumToStringSparse_hh */
